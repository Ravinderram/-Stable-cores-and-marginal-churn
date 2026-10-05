"""Download the monthly river water-quality records used by steps 04, 05b and 07 (sampling variability, fixed
effort, compliance rules).

Source: Esri India Living Atlas, layer "Time Series River Water Quality 2020_2025" (CPCB, CWC and state agency
records, one record per station and month). The records are not redistributed in this repository; this script
recreates data/monthly/esri_timeseries_2020_2025.csv from the public feature service.

The analysis used the copy downloaded on 4 October 2026 (34,037 records; see data/monthly/esri_timeseries_meta.json).
The layer may have changed since, so the script prints the record count and per-agency totals for comparison.
Written for this repository with the standard ArcGIS REST query interface; the build environment could not reach the
service, so check the printed counts against the metadata file before running the pipeline.

Usage (from the repository root):  python code/tools/download_esri_monthly.py
Only the Python standard library is needed.
"""
import csv, datetime, json, os, time, urllib.parse, urllib.request

URL = 'https://livingatlas.esri.in/server1/rest/services/Water/Water_Quality_of_Rivers_in_India/FeatureServer/1'
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
OUT = os.path.join(ROOT, 'data', 'monthly', 'esri_timeseries_2020_2025.csv')


def get(path='', **params):
    params['f'] = 'json'
    url = URL + (f'/{path}' if path else '') + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.load(r)
            if 'error' in d:
                raise RuntimeError(d['error'])
            return d
        except Exception as e:  # network hiccups: retry with back-off
            if attempt == 4:
                raise
            print('retrying after error:', e)
            time.sleep(5 * (attempt + 1))


def main():
    layer = get()
    fields = [f for f in layer['fields'] if f.get('type') != 'esriFieldTypeGeometry']
    names = [f['name'] for f in fields]
    dates = {f['name'] for f in fields if f.get('type') == 'esriFieldTypeDate'}
    oid = layer.get('objectIdField', 'objectid')
    n = get('query', where='1=1', returnCountOnly='true')['count']
    page = min(int(layer.get('maxRecordCount', 1000)), 2000)
    print(f'{layer.get("name")}: {n} records reported by the server; page size {page}')

    rows, offset = [], 0
    while offset < n:
        d = get('query', where='1=1', outFields='*', returnGeometry='false', orderByFields=f'{oid} ASC',
                resultOffset=offset, resultRecordCount=page)
        feats = d.get('features', [])
        if not feats:
            break
        rows += [f['attributes'] for f in feats]
        offset += len(feats)
        print(f'  {offset} / {n}')

    for r in rows:  # date fields come as epoch milliseconds; the analysis reads YYYY-MM from the start of the string
        for k in dates:
            if r.get(k) is not None:
                r[k] = datetime.datetime.fromtimestamp(r[k] / 1000, datetime.timezone.utc).strftime('%Y-%m-%d')

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=names, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    by_agency = {}
    for r in rows:
        by_agency[r.get('agency')] = by_agency.get(r.get('agency'), 0) + 1
    print(f'wrote {len(rows)} records to {OUT}')
    print('records by agency:', json.dumps(by_agency, ensure_ascii=False))
    print('compare with data/monthly/esri_timeseries_meta.json (34,037 records, CPCB 10,001 on 4 October 2026)')


if __name__ == '__main__':
    main()
