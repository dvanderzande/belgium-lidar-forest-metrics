# Belgium LiDAR forest metrics

Forest-structure layers for all of Belgium at 30 m resolution, computed from the regional airborne LiDAR surface and
terrain models (DSM − DTM): canopy cover, 90th percentile height, height standard deviation and mean height. They
reproduce the LiDAR forest covariates of the England red wood ant habitat model
([jwiethase/Ant_ENM_England_2025](https://github.com/jwiethase/Ant_ENM_England_2025)) for use in a Belgian species
distribution model, and come with a quality layer that gives the source survey of every cell and marks a Walloon zone
with unreliable values.

**Status: version 0.1.0.** All of Belgium processed and checked. The output grid (EPSG:31370, cells on multiples of
30 m) is provisional and may change to match the grid of the species distribution model.

- **Download the layers:** the files are attached to the
  [v0.1.0 release](https://github.com/dvanderzande/belgium-lidar-forest-metrics/releases/tag/v0.1.0)
  (GeoTIFFs, zone polygon, data description and checksums; data licence CC BY 4.0 with the source attributions below).
- **Report** (background, methods, decisions, results, data quality, notes for modellers):
  [`docs/Belgium_LiDAR_metrics_report.pdf`](docs/Belgium_LiDAR_metrics_report.pdf)
- **Methods in detail** (definitions, deviations from England, checks, quality layer): [`docs/methods.md`](docs/methods.md)
- **Changes:** [`CHANGELOG.md`](CHANGELOG.md)

## The layers

`metrics_30m_belgium.tif`: 4 bands, float32, NaN outside Belgium, 9400 × 7600 cells of 30 m, EPSG:31370 (Belgian
Lambert 72), tiled and compressed. Every 30 m cell is computed from its 900 canopy-height pixels of 1 m.

| band | name (as in England) | definition |
|---|---|---|
| 1 | `cover_VOM` | share of pixels with data whose height is > 2 m (0–1) |
| 2 | `perc09_height_VOM` | 90th percentile of height (R quantile type 7) |
| 3 | `sd_height_VOM` | standard deviation of height (n − 1) |
| 4 | `mean_height_VOM` | mean height |

The 1 m canopy height is encoded like the Environment Agency Vegetation Object Model used in England: vegetation of at
least 2.5 m keeps its height, everything else (ground, low vegetation, buildings, water) is 0 and counts in all four
metrics. Two deviations from England, needed because the Belgian input is DSM − DTM: canopy gaps up to 3 m are
filled (leaf-off crowns), and heights above 60 m are set to 0 (cranes, wind turbines, masts). See
[`docs/methods.md`](docs/methods.md).

`quality_30m_belgium.tif` (same grid, uint8, NoData 0) and `quality_zone_wallonia_2013_2014.gpkg`:

| code | meaning |
|---|---|
| 1 | DHMV II, Flanders and Brussels (2013–2015, leaf-off) |
| 2 | Wallonia 2013–2014 survey, reliable |
| 3 | Wallonia 2013–2014 survey, zone where leaf-off broadleaf forest is missed: cover and heights too low (3,240 km²) |

## Installation

```bash
conda env create -f environment.yml
conda activate belgium-lidar
cp dev.env.example dev.env      # then set BELGIUM_LIDAR_WORK_DIR in dev.env
source dev.env
python -m unittest discover tests
```

The scripts run from the repository folder; no package installation is needed. A working folder with about 170 GB of
free space is needed for Belgium (280 GB while the Walloon ZIPs are still there).

## Producing the layers

```bash
source dev.env
W="$BELGIUM_LIDAR_WORK_DIR"

# 1. Downloads (server-bound; resumable, finished files are skipped)
python scripts/download_flanders.py "$W/inputs/flanders_dhmv2"          # DHMV II DSM/DTM, GRB buildings: 56 GB, ~6 h
python scripts/download_brussels.py "$W/inputs/brussels_dhmv2"          # DHMV II, UrbIS buildings: 1 GB, ~7 min
python scripts/download_wallonia.py "$W/inputs/wallonia_2013_2014"      # MNS/MNT 2013-14, PICC: 113 GB, ~11 h

# 2. Canopy height and the four metrics, per 6 km tile (resumable), then one merged GeoTIFF
python scripts/process_region.py "$W/inputs" "$W/belgium"               # ~6 h in one process (Walloon tiles ~36 s each)
#    or side by side: run "--part K 6" for K = 1..6 in six terminals, then once without --part to merge (~55 min)

# 3. Checks: grid, value ranges, seams, map; steps at the regional borders
python scripts/check_metrics.py "$W/belgium/metrics_30m_belgium.tif"
python scripts/check_borders.py "$W/belgium/metrics_30m_belgium.tif" "$W/inputs"

# 4. Quality layer (needs the SPW forest mask 2018)
python scripts/download_forest_mask.py "$W/inputs"                      # 10 MB, ~2 min
python scripts/make_quality_mask.py "$W/belgium/metrics_30m_belgium.tif" "$W/inputs"
```

One region only: `--region flanders|brussels|wallonia` (repeatable); one test tile: `--tile X Y` (lower-left corner).
Every output GeoTIFF records its sources, rules and script in its tags.

**Evidence behind the Walloon zone** (optional; reproduces the analysis in [`docs/methods.md`](docs/methods.md)):

```bash
python scripts/download_border_strip.py "$W/inputs" "$W/inputs/dhmv2_wallonia_strip"   # 3.3 GB, ~18 min
python scripts/compare_border_strip.py "$W/inputs" "$W/border_strip"
python scripts/check_wallonia_quality.py "$W/belgium/metrics_30m_belgium.tif" "$W/inputs" "$W/wallonia_quality"
```

## Code

| part | what it does |
|---|---|
| `belgium_lidar_metrics/canopy.py` | canopy height rules (2.5 m, 60 m, buildings, 3 m gap filling) |
| `belgium_lidar_metrics/metrics.py` | the four England metrics per 30 m cell |
| `belgium_lidar_metrics/download.py`, `output.py` | WCS/WFS downloads with retries and a size limit; GeoTIFF writing with provenance tags |
| `scripts/download_*.py` | input downloads per region, the SPW forest mask and the border strip |
| `scripts/process_region.py` | canopy height and metrics per tile for any regions, then the merged map |
| `scripts/check_*.py`, `compare_border_strip.py` | checks of the result and of the Walloon data quality |
| `scripts/make_quality_mask.py` | the quality layer and the zone polygon |
| `tests/` | unit tests of the rules and metrics against hand-computed values (incl. NoData and all-zero cells) |

## Data sources and attribution

The repository contains no data; the scripts download it from the public services of the three regions.

| data | source | licence / attribution |
|---|---|---|
| DHMV II DSM and DTM 1 m, GRB buildings, Flemish border | Digitaal Vlaanderen | Gratis open data licentie Vlaanderen v1.0; « Digitaal Vlaanderen » |
| UrbIS buildings and Brussels border | Brussels Region (UrbIS, Paradigm) | open data; CC-0 stated in the catalogue (to confirm for the buildings) |
| MNS and MNT 2013–2014, PICC buildings, Walloon border, forest mask 2018 | Service public de Wallonie | free reuse; « Sources des données : SPW » on every map or document |

Maps and documents made from the layers must carry these attributions.

## Relation to the England model

The metric definitions follow the R code of the England model (Joris Wiethase, `Ant_ENM_England_2025`), which is the
reference for this project. This repository is an independent Python implementation: it contains no code from that
repository. Please cite the England model as well when you use these layers in the same way.

## Citation, licence and acknowledgement

- Cite this repository with the information in [`CITATION.cff`](CITATION.cff) (GitHub: "Cite this repository").
- Code under the [MIT licence](LICENSE); the released layers under CC BY 4.0, with the source attributions above.
- Parts of the code and documentation were written with the help of an AI coding assistant. All methods, results and
  checks were reviewed by the author.
- Questions and problems: please open an issue on GitHub.

Author: Dimitry Van der Zande (Institute of Natural Sciences).
