"""Rev 20 local step 07: river size and upstream population for located monitoring stations (reviewer Major Comment 5).

Run on your own computer (the cloud workspace cannot download HydroSHEDS or WorldPop). Same virtual environment as the
earlier local kit (geopandas, rasterio, numpy, pandas, requests).

What it does
  1. Downloads HydroRIVERS v1.0 Asia and HydroBASINS v1c Asia level 12 (or reuses copies already on disk).
  2. Snaps each station to a HydroRIVERS reach: among reaches within 1 km the one with the largest upstream area,
     otherwise the nearest reach within 2 km; stations farther than 2 km are left unsnapped (never guessed).
  3. Long-term average discharge (DIS_AV_CMS, modelled, 1971-2000) and upstream area of the snapped reach.
  4. Population: WorldPop 100 m constrained rasters you already downloaded for the population step (India, 2016-2024)
     are summed to 1-km cells, cells are assigned to HydroBASINS level-12 units, and population is accumulated
       (a) over the whole upstream catchment (all level-12 units draining to the station's unit, own unit included), and
       (b) over level-12 units drained by reaches within 50 km upstream along the river network (own unit included).
     Population outside India is not counted (India rasters only); transboundary catchments are flagged.
  5. Population per unit discharge (people per m3/s) as a dilution proxy.
Output: data/context/upstream_station.csv (one row per station) and data/context/upstream_report.json.
Nothing in the water-quality data is changed.

Usage (from the repository root):
    python code/local/07_upstream_hydrorivers.py
Downloads go to code/local/hydrosheds/ (git-ignored). Put the WorldPop rasters in code/local/population/ or edit
WORLDPOP_DIR below.
"""
import json, math, sys, time, zipfile
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
KIT = HERE.parent
OUT = KIT.parent / 'data' / 'context'; OUT.mkdir(parents=True, exist_ok=True)
DATA = HERE / 'hydrosheds'; DATA.mkdir(exist_ok=True)
STATIONS = HERE / 'stations_located_r20.csv'
# folder holding ind_pop_<year>_CN_100m_R2025A_v1.tif from the population step; searched recursively if not found
WORLDPOP_DIR = HERE / 'population'
POP_YEAR = 2020                       # one mid-period year; the series is smooth (median growth 6.9% over 2016-2024)
BBOX = (60.0, 5.0, 105.0, 40.0)       # lon/lat window that holds all Indian catchments incl. upstream parts abroad
URLS = {'rivers': 'https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_as_shp.zip',
        'basins': 'https://data.hydrosheds.org/file/hydrobasins/standard/hybas_as_lev12_v1c.zip'}

def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)

def download(url, dest):
    import requests
    if dest.exists() and dest.stat().st_size > 1e6:
        log('present:', dest.name); return dest
    tmp = dest.with_suffix(dest.suffix + '.part'); pos = tmp.stat().st_size if tmp.exists() else 0
    for attempt in range(8):
        try:
            h = {'Range': f'bytes={pos}-'} if pos else {}
            with requests.get(url, headers=h, stream=True, timeout=120) as r:
                if r.status_code == 416: break
                r.raise_for_status()
                mode = 'ab' if pos and r.status_code == 206 else 'wb'
                if mode == 'wb': pos = 0
                with open(tmp, mode) as f:
                    for ch in r.iter_content(1 << 20):
                        f.write(ch); pos += len(ch)
            break
        except Exception as e:
            log('download interrupted, retrying:', e); time.sleep(10 * (attempt + 1))
    tmp.rename(dest); log('downloaded', dest.name, round(dest.stat().st_size / 1e6), 'MB'); return dest

def find_shp(pattern):
    hits = sorted(KIT.rglob(pattern))
    return hits[0] if hits else None

def get_layer(kind, pattern):
    shp = find_shp(pattern)
    if shp: log('using existing', shp); return shp
    z = download(URLS[kind], DATA / Path(URLS[kind]).name)
    with zipfile.ZipFile(z) as zz: zz.extractall(DATA)
    shp = find_shp(pattern)
    if not shp: sys.exit(f'could not find {pattern} after unzipping {z}')
    return shp

def main():
    import geopandas as gpd, rasterio
    from rasterio.features import rasterize
    from rasterio.enums import Resampling
    S = pd.read_csv(STATIONS)
    log('stations:', len(S))
    riv_shp = get_layer('rivers', 'HydroRIVERS_v10_as.shp')
    bas_shp = get_layer('basins', 'hybas_as_lev12_v1c.shp')
    log('reading rivers (bbox)'); R = gpd.read_file(riv_shp, bbox=BBOX)
    R = R[['HYRIV_ID', 'NEXT_DOWN', 'MAIN_RIV', 'LENGTH_KM', 'DIST_DN_KM', 'UPLAND_SKM', 'DIS_AV_CMS', 'ORD_STRA', 'HYBAS_L12', 'geometry']]
    log('reaches:', len(R))
    log('reading level-12 basins (bbox)'); Bz = gpd.read_file(bas_shp, bbox=BBOX)[['HYBAS_ID', 'NEXT_DOWN', 'UP_AREA', 'SUB_AREA', 'geometry']]
    log('level-12 units:', len(Bz))

    # ---- 1. snap stations to reaches (metric CRS for distances)
    P = gpd.GeoDataFrame(S, geometry=gpd.points_from_xy(S.lon, S.lat), crs=4326).to_crs(7755)   # WGS 84 / India NSF LCC
    Rm = R.to_crs(7755)
    sidx = Rm.sindex
    rows = []
    for i, p in P.iterrows():
        cand = list(sidx.query(p.geometry.buffer(2000)))
        rec = dict(station_uid=p.station_uid, snap_rule='unsnapped (> 2 km from any reach)')
        if cand:
            c = Rm.iloc[cand].copy(); c['d'] = c.distance(p.geometry)
            c = c[c.d <= 2000]
            if len(c):
                near = c[c.d <= 1000]
                pick = near.sort_values('UPLAND_SKM', ascending=False).iloc[0] if len(near) else c.sort_values('d').iloc[0]
                rule = 'largest upstream area within 1 km' if len(near) else 'nearest reach within 2 km'
                rec.update(snap_rule=rule, HYRIV_ID=int(pick.HYRIV_ID), snap_dist_m=round(float(pick.d)), UPLAND_SKM=float(pick.UPLAND_SKM),
                           DIS_AV_CMS=float(pick.DIS_AV_CMS), ORD_STRA=int(pick.ORD_STRA), HYBAS_L12=int(pick.HYBAS_L12),
                           MAIN_RIV=int(pick.MAIN_RIV), DIST_DN_KM=float(pick.DIST_DN_KM), n_reaches_within_1km=int(len(near)))
        rows.append(rec)
    SN = pd.DataFrame(rows)
    log('snapped:', int(SN.HYRIV_ID.notna().sum()), 'of', len(SN))

    # ---- 2. population per level-12 unit on a 1-km grid
    tifs = sorted(WORLDPOP_DIR.rglob(f'ind_pop_{POP_YEAR}_CN_100m_R2025A_v1.tif')) if WORLDPOP_DIR.exists() else []
    if not tifs: tifs = sorted(KIT.rglob(f'ind_pop_{POP_YEAR}_CN_100m_R2025A_v1.tif'))
    if not tifs: sys.exit(f'WorldPop raster ind_pop_{POP_YEAR}_CN_100m_R2025A_v1.tif not found under {KIT}; set WORLDPOP_DIR')
    tif = tifs[0]; log('population raster:', tif)
    with rasterio.open(tif) as src:
        f = 10                                                   # 10 x 10 cells of 3 arc-s = 30 arc-s (~1 km)
        h, w = src.height // f, src.width // f
        pop1k = np.zeros((h, w), dtype='float64')
        for r0 in range(0, h, 50):                               # block-wise to keep memory low
            rr = min(50, h - r0)
            win = rasterio.windows.Window(0, r0 * f, w * f, rr * f)
            a = src.read(1, window=win, masked=True).filled(0).astype('float64')
            a[a < 0] = 0
            pop1k[r0:r0 + rr] = a.reshape(rr, f, w, f).sum(axis=(1, 3))
        tr = src.transform * src.transform.scale(f, f)
        total_pop = pop1k.sum()
    log('population on 1-km grid:', round(total_pop / 1e6, 1), 'million')
    Bi = Bz.to_crs(4326).reset_index(drop=True)
    Bi['idx'] = np.arange(1, len(Bi) + 1)
    lab = rasterize(((g, v) for g, v in zip(Bi.geometry, Bi.idx)), out_shape=pop1k.shape, transform=tr, fill=0, dtype='int32')
    unit_pop = np.bincount(lab.ravel(), weights=pop1k.ravel(), minlength=len(Bi) + 1)
    Bi['pop'] = unit_pop[Bi.idx.values]
    log('population assigned to level-12 units:', round(Bi['pop'].sum() / 1e6, 1), 'million (', round(100 * Bi['pop'].sum() / total_pop, 1), '% )')
    popmap = dict(zip(Bi.HYBAS_ID.astype('int64'), Bi['pop']))

    # ---- 3. upstream accumulation on the level-12 graph (NEXT_DOWN)
    up = {}
    for hid, nd in zip(Bi.HYBAS_ID.astype('int64'), Bi.NEXT_DOWN.astype('int64')):
        up.setdefault(nd, []).append(hid)
    def upstream_units(start):
        seen, stack = {start}, [start]
        while stack:
            for u in up.get(stack.pop(), []):
                if u not in seen: seen.add(u); stack.append(u)
        return seen
    # reach graph for the 50-km window
    rup = {}
    for rid, nd in zip(R.HYRIV_ID.astype('int64'), R.NEXT_DOWN.astype('int64')):
        rup.setdefault(nd, []).append(rid)
    rlen = dict(zip(R.HYRIV_ID.astype('int64'), R.LENGTH_KM)); rbas = dict(zip(R.HYRIV_ID.astype('int64'), R.HYBAS_L12.astype('int64')))
    def units_within(rid, km=50.0):
        units, stack = {rbas[rid]}, [(rid, 0.0)]
        while stack:
            r, dist = stack.pop()
            for u in rup.get(r, []):
                dd = dist + rlen.get(u, 0)
                if dd <= km: units.add(rbas[u]); stack.append((u, dd))
        return units
    outside = set(Bi.HYBAS_ID[Bi['pop'] == 0].astype('int64'))
    res = []
    for _, r in SN.iterrows():
        rec = r.to_dict()
        if pd.notna(r.get('HYRIV_ID')):
            u12 = int(r.HYBAS_L12)
            U = upstream_units(u12); U50 = units_within(int(r.HYRIV_ID))
            rec.update(pop_upstream=round(sum(popmap.get(u, 0) for u in U)), n_units_upstream=len(U),
                       pop_upstream_50km=round(sum(popmap.get(u, 0) for u in U50)), n_units_50km=len(U50),
                       pop_own_unit=round(popmap.get(u12, 0)),
                       share_units_without_india_pop=round(len(U & outside) / len(U), 3))
            rec['pop_per_cms'] = round(rec['pop_upstream'] / r.DIS_AV_CMS, 1) if r.DIS_AV_CMS and r.DIS_AV_CMS > 0 else None
        res.append(rec)
    OUTF = pd.DataFrame(res)
    OUTF['pop_year'] = POP_YEAR
    OUTF['source'] = 'HydroRIVERS v1.0 / HydroBASINS v1c level 12 (Lehner & Grill 2013); WorldPop Global2 R2025A constrained 100 m, India only'
    OUTF.to_csv(OUT / 'upstream_station.csv', index=False)
    rep = dict(stations=len(S), snapped=int(OUTF.HYRIV_ID.notna().sum()), snap_rules=OUTF.snap_rule.value_counts().to_dict(),
               median_snap_dist_m=float(OUTF.snap_dist_m.median()), pop_year=POP_YEAR, raster=str(tif.name),
               population_grid_total_million=round(total_pop / 1e6, 1),
               population_in_level12_units_million=round(Bi['pop'].sum() / 1e6, 1),
               stations_with_any_upstream_unit_outside_india=int((OUTF.share_units_without_india_pop > 0).sum()))
    (OUT / 'upstream_report.json').write_text(json.dumps(rep, indent=1))
    log('wrote', OUT / 'upstream_station.csv'); print(json.dumps(rep, indent=1))

if __name__ == '__main__':
    main()
