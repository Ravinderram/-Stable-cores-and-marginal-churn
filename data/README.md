# Data

All files are UTF-8 CSV unless stated. Original reported values are never modified; every parsed value carries a
status flag, and every exclusion or merge is recorded row by row.

## raw/

| File | Content |
|---|---|
| `water_quality_merged_2016-2024.xlsx` | The tables of the nine CPCB annual river water-quality reports (2016-2024) in one workbook: 11,157 station-year rows. Sheets: `Merged Data` (values as printed, with source file, sheet and row), `Summary`, `Processing Log`. Never edited. |

## processed/

| File | Content |
|---|---|
| `integrated_station_year_enriched.csv` | One row per workbook row (11,157; 120 columns). Main input of the pipeline. |
| `station_dictionary_r20.csv` | One row per harmonized station (1,891): codes, original names, state, years observed, harmonization confidence, exclusion reasons, coordinates, 2023 sampling frequency. |
| `row_to_station_map_r20.csv` | Workbook row to harmonized station, with source file and source row. |

Column groups of `integrated_station_year_enriched.csv`:

| Columns | Meaning |
|---|---|
| `row_id`, `Year`, `Station Code`, `Monitoring Location`, `State`, `Table Title`, `Source File`, `Source Row` | Identification of the workbook row |
| `Temperature (°C) Min` ... `Fecal Streptococci (MPN/100 mL) Max` | Values exactly as printed in the reports |
| `TEMP_min`, `DO_min`, `PH_max`, `COND_max`, `BOD_max`, `NIT_max`, `FC_max`, `TC_max`, `FS_max` (min and max for each) | Parsed numeric values |
| `<value>_status` | `numeric`, `BDL` (below detection limit), `missing` or `qc_excluded` |
| `station_uid`, `harmonized_name`, `harmonization_confidence` | Harmonized station identity; confidence `high`, `medium`, `low`, `ambiguous` or `single_year` |
| `use_in_longitudinal`, `exclusion_reason` | Longitudinal-use flag; reasons: single-year station, duplicate row with identical values, station listed twice in a year with different values |
| `lat`, `lon`, `coord_confidence`, `coord_source` | Coordinates from official listings only (never guessed); confidence `High`, `Medium`, `Low` or `Ambiguous / unresolved` |
| `reported_frequency`, `frequency_reference` | Sampling frequency reported in the CPCB 2023 table (monthly, quarterly, yearly); a station attribute that refers to 2023 |
| `river_cpcb_table`, `basin_cpcb_table`, `river_from_name`, `district` | River, basin and district as given or parsed from the name |
| `rain_*`, `imd_cell_*`, `heavy_rain_days`, ... | IMD rainfall of the nearest 0.25 degree cell (located stations) |
| `pop_*` | WorldPop population within 1 km and 5 km (located stations) |
| `hybas_*`, `basin__*` | HydroBASINS level 4 and level 6 assignment |
| `state_*_2021inv`, `state_stp_reference` | State sewage generation and treatment capacity (CPCB National Inventory of STPs, March 2021) |
| `season_of_min_max` | Always "Season unavailable": the reports do not say when the minimum or maximum was sampled |

## context/

| File | Content |
|---|---|
| `rainfall_station_year.csv` | IMD rainfall per located station and year: annual and seasonal totals, heavy-rain days, anomalies against 1991-2020. |
| `population_station_year.csv` | WorldPop population and density within 1 km and 5 km per located station and year. |
| `basin_station.csv`, `basin_disagreements_for_review.csv` | HydroBASINS level 4 and 6 assignments (point in polygon), and the stations whose main basin differs from the usual one for their CPCB basin name. |
| `upstream_station.csv`, `upstream_report.json` | HydroRIVERS snapping (1,171 of 1,187 located stations within 2 km; median 307 m), modelled long-term discharge, upstream area, upstream-catchment population and population per unit discharge. Produced by `code/local/07_upstream_hydrorivers.py`. India-only population: 348 stations have upstream units without Indian population. |

## monthly/

| File | Content |
|---|---|
| `esri_monthly_links.csv` | One-to-one links between Esri India monitoring keys (`agency|station|lat|lon`) and harmonized stations (290 linked). |
| `esri_timeseries_meta.json` | Description of the monthly layer as downloaded on 4 October 2026 (fields, record counts by agency and year). |
| `esri_timeseries_2020_2025.csv` | **Not included.** Create it with `python code/tools/download_esri_monthly.py`. Only the CPCB records of 2020 and 2021 are used. |

## geo/

| File | Content |
|---|---|
| `india_states_ne10m.gpkg`, `neighbours_ne10m.gpkg` | Natural Earth 1:10m polygons of India's states (admin-1) and of neighbouring countries (admin-0), used only to draw Fig. 2. |
| `station_use_r20.csv` | Optional. If present, step 08 takes the main spatial sample from it; otherwise the sample is the stations with High or Medium coordinate confidence that lie within 2 km of a HydroRIVERS reach (1,155 stations). |

## analysis/ (built by the pipeline)

| File | Built by | Content |
|---|---|---|
| `analysis_r20.csv` | steps 00 and 02 | The enriched table plus panel flags (`in_panel_A`, `in_panel_A_long`), cleaned state, adverse status (`ADV_DO`, `ADV_BOD`, `ADV_PH`, `ADV_FC`, `ADV_FC_cens`, `ADV_MULTI_GE2`, `ADV_MULTI_GE2_cens`, `ADV_ORG` = BOD or DO adverse) |
| `stations_r20.csv` | step 08 | Station covariates for the context models, including `spatial_main` |
| `context_r20.csv` | step 08 | Station-year table for the context models |

`ADV_*` is 1 (adverse), 0 (compliant) or missing. `ADV_FC` codes FC as in the earlier draft (ceiling states and
1,600 compliant); `ADV_FC_cens` is the primary coding, with FC missing where compliance cannot be assessed.
