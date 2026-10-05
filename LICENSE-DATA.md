# Data licence

The harmonized datasets, station dictionary, links, derived station-level variables and result tables produced in
this project (`data/processed/`, `data/context/`, `data/monthly/esri_monthly_links.csv`, `analysis/`, `results/`,
`figures/`) are released under the **Creative Commons Attribution 4.0 International licence (CC BY 4.0)**,
https://creativecommons.org/licenses/by/4.0/. Please cite the data deposit and the article (see `CITATION.cff`).

## Third-party sources

The derived variables keep the attribution requirements of their sources. Cite them when you reuse those variables.

| Source | Terms | Used in |
|---|---|---|
| Central Pollution Control Board, annual river water-quality reports 2016-2024 | Public reports of the Government of India; values are reproduced as printed, with source file and row | `data/raw/`, `data/processed/` |
| Esri India Living Atlas, "Time Series River Water Quality 2020_2025" | Not redistributed here; download from the public feature service under its own terms | `data/monthly/` (links only) |
| India Meteorological Department, gridded daily rainfall (0.25 degree; Pai et al., 2014, Mausam 65(1)) | Cite IMD and Pai et al. (2014) | `data/context/rainfall_station_year.csv` |
| WorldPop, Global2 constrained population, 100 m, release 2025A | CC BY 4.0; cite WorldPop (University of Southampton) | `data/context/population_station_year.csv`, `data/context/upstream_station.csv` |
| HydroSHEDS: HydroRIVERS v1.0 and HydroBASINS v1c | See the HydroSHEDS licence agreement at https://www.hydrosheds.org; cite Lehner and Grill (2013) | `data/context/basin_station.csv`, `data/context/upstream_station.csv` |
| Natural Earth, admin-0 and admin-1, 1:10m | Public domain | `data/geo/` |

Map lines in the figures delineate study areas and do not necessarily depict accepted national boundaries.
