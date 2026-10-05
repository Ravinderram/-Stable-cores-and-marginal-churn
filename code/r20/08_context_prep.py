"""Rev 20, step 8a: station-year table for the context models (rainfall, population, basin, sampling frequency,
selection into the located sample).
Main spatial sample: stations with use_recommendation == 'main' in Phase G (data/geo/station_use_r20.csv) when that file is
present; otherwise High or Medium coordinate confidence (fallback, flagged in the output).
Rainfall: annual IMD total of the nearest 0.25 deg cell; exposure = ln(year total / cell mean 2016-2024); climate =
ln(cell mean 2016-2024). Flagged cells: |z| > 2 of the cell's mean 2016-2024 anomaly against 1991-2020.
Population: log10 of the station's mean 2016-2024 WorldPop density within 5 km (1 km alternative).
Output: analysis/context_r20.csv, analysis/stations_r20.csv"""
import sys, os
from core import *
d = load(None)
R = pd.read_csv(f'{B}/data/context/rainfall_station_year.csv')
R['cell'] = R.imd_cell_lat.round(3).astype(str) + '_' + R.imd_cell_lon.round(3).astype(str)
cm = R.groupby('cell').rain_annual_mm.mean(); ca = R.groupby('cell').rain_annual_anomaly_pct.mean()
z = (ca - ca.mean()) / ca.std(); flagged = set(z.index[z.abs() > 2])
R['rain_clim'] = np.log(R.cell.map(cm)); R['rain_expo'] = np.log(R.rain_annual_mm / R.cell.map(cm)); R['rain_flag'] = R.cell.isin(flagged)
P = pd.read_csv(f'{B}/data/context/population_station_year.csv')
pm = P.groupby('station_uid')[['pop_density_5km_per_km2', 'pop_density_1km_per_km2']].mean()
use_file = f'{B}/data/data/geo/station_use_r20.csv'
st = d[d.in_panel_A].groupby('station_uid').agg(state=('station_state', 'first'), coord=('coord_confidence', 'first'),
                                                lat=('lat', 'first'), lon=('lon', 'first'), basin_sub=('hybas_lev06', 'first'),
                                                basin_main=('hybas_main_basin_lev04', 'first'), n_obs=('Year', 'nunique'), first=('Year', 'min'),
                                                freq=('reported_frequency', lambda s: s.dropna().mode().iat[0] if s.notna().any() else 'not reported'),
                                                bod_med=('BOD_max', 'median'), cond_med=('COND_max', 'median'), temp_med=('TEMP_max', 'median'))
if os.path.exists(use_file):
    U = pd.read_csv(use_file).set_index('station_uid'); st['spatial_main'] = U.reindex(st.index).spatial_main.fillna(False).astype(bool)
    st['spatial_rule'] = 'Phase G use_recommendation'
else:
    # rev 20 rule, reproducible from deposited fields: High or Medium coordinate confidence and within 2 km of a
    # HydroRIVERS reach (local step 07); replaces the Phase G 'all checks passed' flag, which is not in the deposit
    snapped = set(pd.read_csv(f'{B}/data/context/upstream_station.csv').query('HYRIV_ID == HYRIV_ID').station_uid) \
        if os.path.exists(f'{B}/data/context/upstream_station.csv') else set(st.index)
    st['spatial_main'] = st.coord.isin(['High', 'Medium']) & st.index.isin(snapped)
    st['spatial_rule'] = 'High or Medium confidence and within 2 km of a HydroRIVERS reach'
st['located'] = st.lat.notna()
st = st.join(pm)
st['lpop5'] = np.log10(st.pop_density_5km_per_km2.clip(lower=0.5)); st['lpop1'] = np.log10(st.pop_density_1km_per_km2.clip(lower=0.5))
up_file = f'{B}/data/context/upstream_station.csv'
if os.path.exists(up_file):
    UP = pd.read_csv(up_file).set_index('station_uid')
    st = st.join(UP[['HYRIV_ID', 'snap_dist_m', 'UPLAND_SKM', 'DIS_AV_CMS', 'pop_upstream', 'pop_upstream_50km', 'pop_per_cms', 'share_units_without_india_pop']])
st.to_csv(f'{B}/analysis/stations_r20.csv')
C = d[d.in_panel_A_long].merge(R[['station_uid', 'year', 'cell', 'rain_clim', 'rain_expo', 'rain_flag']].rename(columns={'year': 'Year'}),
                                on=['station_uid', 'Year'], how='left')
keep = ['station_uid', 'Year', 'ADV_BOD', 'ADV_DO', 'ADV_PH', 'ADV_FC', 'ADV_FC_cens', 'ADV_MULTI_GE2', 'ADV_MULTI_GE2_cens', 'cell', 'rain_clim', 'rain_expo', 'rain_flag']
C = C[keep].merge(st[['state', 'basin_sub', 'basin_main', 'spatial_main', 'located', 'lpop5', 'lpop1', 'freq', 'n_obs', 'first'] +
                     ([c for c in ['DIS_AV_CMS', 'pop_upstream', 'pop_upstream_50km', 'pop_per_cms'] if c in st.columns])], left_on='station_uid', right_index=True)
C.to_csv(f'{B}/analysis/context_r20.csv', index=False)
print(st.spatial_rule.iat[0], '| spatial main stations', int(st.spatial_main.sum()), 'station-years (A_long)', int(C.spatial_main.sum()),
      '| flagged cells', len(flagged), '| SD rain_expo', round(C.rain_expo[C.spatial_main].std(), 3), '-> x', round(np.exp(C.rain_expo[C.spatial_main].std()), 3),
      '| SD rain_clim (stations)', round(st.join(C.groupby('station_uid').rain_clim.first()).rain_clim[st.spatial_main].std(), 3))
