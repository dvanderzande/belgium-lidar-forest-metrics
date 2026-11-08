# Methods

How the Belgian 30 m forest-structure layers are made, why they are made this way, and how they were checked. The
report (`docs/Belgium_LiDAR_metrics_report.pdf`) gives the same content with figures.

## 1. Reference: the England wood-ant model

The layers reproduce the LiDAR forest covariates of the England red wood ant habitat model (Joris Wiethase,
[`Ant_ENM_England_2025`](https://github.com/jwiethase/Ant_ENM_England_2025), reviewed at commit 4792f5b). Paths below
are in that repository.

**Input in England:** the Environment Agency **Vegetation Object Model (VOM)**, a 1 m raster of vegetation height
above ground with buildings removed, derived from the national LiDAR programme. Only objects classified as vegetation
and at least 2.5 m tall are included.

**Metrics** per 30 m cell from its ~900 one-metre pixels (`R/HPC_VOM_scripts/gridAggregate_VOM_HPC_V2.R`, functions
in `source/misc_functions.R`), all ignoring NoData (`na.rm = TRUE`):

| layer | definition | use in the England model |
|---|---|---|
| `cover_VOM` | share of non-NoData pixels with height > 2 m (0–1) | forest PCA, forest mask, forest patches |
| `perc09_height_VOM` | 90th percentile of height (R `quantile`, type 7) | forest PCA |
| `sd_height_VOM` | standard deviation of height (n − 1) | forest PCA |
| `mean_height_VOM` | mean height | computed, not used |

**Use downstream:**
- `R/S4b_process_forestTopoEffort_local.R`: forest = cover ≥ 0.3; forest buffered by 70 m is the model domain (95 % of
  nests lie within 67 m of forest). Distance to forest was computed but not kept as a covariate.
- `R/S5_run_forest_PCA_HPC.R`: PCA (centred, scaled) of cover, √p90 and √SD inside the domain gives `forest_PC1` and
  `forest_PC2`. Only the loadings are saved, not the centring and scaling.
- `R/S11_get_patch_distance_HPC.R`: cover ≥ 0.3 again for forest patches and their area (suitable, colonisable and
  translocation sites).
- The rasters were reprojected with `terra::project` (bilinear by default).

Because cover enters as a threshold and as patch areas, and p90 and SD through a PCA, the Belgian layers must be
encoded exactly as in England.

### No-vegetation encoding in the EA VOM

Checked on two 5 km tiles of the 2018 VOM (Forest of Dean, SO51se and SO61sw; flown 30 Nov 2017 – 20 Apr 2018,
leaf-off): there are **no NoData pixels**; ground, fields, roads and buildings are **0**, and every other pixel is
**≥ 2.5 m** (nothing between 0 and 2.5 m). The England code aggregates the tiles as they are, so the zeros count in
all four metrics. With NoData instead of zeros, cover would be 1 in every cell with any vegetation.

| 30 m cells (England definitions) | SO51se (forest and farmland) | SO61sw (mostly forest) |
|---|---|---|
| cover, median (5–95 %) | 0.27 (0–0.97) | 0.81 (0.01–0.99) |
| p90, median (5–95 %) | 7.6 m (0–28.6) | 20.1 m (0–31.2) |
| forest cells (cover ≥ 0.3) | 48.9 % | 85.2 % |
| in forest cells, median p90 / SD / mean | 20.6 / 6.7 / 11.3 m | 21.1 / 6.9 / 12.3 m |

These serve as a plausibility reference, not a calibration target (a different forest).

## 2. Belgian input data

| region | heights (1 m) | buildings | border |
|---|---|---|---|
| Flanders | DHMV II DSM and DTM (Digitaal Vlaanderen), flown 2013–2015 leaf-off (e.g. Sonian Forest: 25 Mar – 1 Apr 2014, ~5.7 points/m²); WCS in 1 km pieces | GRB buildings (GBG) and annexes (GBA) | VRBG `Refgew` |
| Brussels | DHMV II (the Flemish survey covers Brussels completely) | UrbIS Buildings (249,072) | UrbIS `Region` |
| Wallonia | SPW MNS and MNT 2013–2014, flown 12/2012 – 03/2014, ~1.5 points/m²; one GeoTIFF per province | PICC building footprints (3.87 M) | SPW administrative limits, layer 0 |
| check only | SPW forest mask 2018 (Gembloux Agro-Bio Tech / ULiège; woody vegetation ≥ ~2 m from orthophotos 2006–2018) | | |

All in EPSG:31370 (the forest mask is reprojected from EPSG:3812). Both height sources are leaf-off raster DSM/DTM
products of the same period and resolution, which keeps the logic of the England VOM.

**Not used:** the LAZ point clouds (about 1.5–3 TB for Flanders alone), and the Walloon survey of 2021–2022 (0.5 m,
~13 points/m², about 400 GB, eight years later than DHMV II).

**Regional borders:** the three borders fit within centimetres (overlap Flanders–Wallonia 1.45 ha, hairline gaps of
7 cm on average). Holes in their union are real: German enclaves along the Vennbahn and Dutch enclaves in
Baarle-Hertog.

## 3. Processing

**Canopy height at 1 m** (`belgium_lidar_metrics/canopy.py`), in this order:
1. height = DSM − DTM;
2. heights above 60 m → 0 (not vegetation; first, so that gap filling cannot spread them);
3. building footprints → 0;
4. gap filling: grey-scale closing with a 3 × 3 m window (a maximum, then a minimum filter) fills holes up to about
   2 m wide in crowns without widening crowns or forest edges; buildings are set to 0 again afterwards;
5. heights below 2.5 m → 0 (the EA VOM encoding);
6. NaN only where there are no height or building data.

**Metrics per 30 m cell** (`belgium_lidar_metrics/metrics.py`): the four England definitions above; zeros count, NaN
is ignored; SD is NaN in a cell with fewer than 2 pixels with data (as R).

**Grid and tiling** (`scripts/process_region.py`): 30 m cells on multiples of 30 m in EPSG:31370 (provisional; one
setting, `GRID_ORIGIN`). The work goes per 6 km tile, read with a 10 m halo so that gap filling and cells across
download blocks come out as in one large raster. Where a tile crosses a regional border, each region's canopy height
is made from its own heights and buildings inside its own border, and the 1 m heights are joined before the 30 m
cells are computed, so cells on a regional border are complete. A pixel inside two borders goes to the region listed
first (Flanders, Brussels, Wallonia).

## 4. Deviations from England

### 4.0 The input: an approximation of the EA Vegetation Object Model

England computed the metrics from the EA VOM, a 1 m raster that the Environment Agency derived from the classified
point clouds of its national LiDAR programme: only points classified as vegetation, buildings removed, vegetation
modelled into discrete objects of at least 2.5 m (crowns filled). Belgium has no such product, so the Belgian input is
an approximation built from the regional surface and terrain models. The metric definitions are identical.

| | England (EA VOM) | Belgium (this work) |
|---|---|---|
| source | classified point clouds → VOM raster (by the EA) | DSM and DTM rasters of the regions, 1 m |
| what counts as vegetation | points classified as vegetation | every height except buildings (footprint maps) and heights above 60 m |
| tree crowns | modelled into whole objects | gaps up to 3 m filled (grey-scale closing, 4.1) |
| below 2.5 m | 0 | 0 (same rule) |
| metrics per 30 m cell | England R code | same definitions (Python, unit-tested) |

**Why the VOM approach was not applied:** reproducing the VOM would mean reprocessing the raw point clouds
(classifying vegetation points and modelling tree objects). That is too big a job for this delivery: about 1.5–3 TB
of point clouds for Flanders alone (and about 850 GB for the Walloon survey of 2021–2022), and the EA's object
modelling is not published in enough detail to be reproduced exactly.

**Effect of the approximation:** cover and the 0.3 forest threshold are comparable with England; heights in forest
are 3–6 m higher (gap filling, 4.1); structures that are not buildings and are lower than 60 m (pylons, overhead
lines, cranes, industrial installations) remain in a few places, mainly ports and industrial zones (4.2). A mask from
an independent vegetation map could remove most of them in a later version. The band names keep the England names
(`…_VOM`) so that the England scripts can be reused unchanged, although the layers are not computed from the EA VOM.

### 4.1 Gap filling of leaf-off crowns (3 m closing)

In the leaf-off DHMV II survey the laser reaches the ground through bare deciduous crowns. In a pilot in the Sonian
Forest (beech, flown late March 2014), 40 % of the pixels inside tall forest were 0 and cover was 0.60, against 14 %
and 0.90 in the EA VOM. The unfilled 1 m map is speckled with ground hits in deciduous stands, while compact (conifer)
stands are solid. A grey-scale closing before the 2.5 m rule was tested:

| Sonian Forest pilot | no filling | closing 3 m | closing 5 m | Forest of Dean (EA VOM) |
|---|---|---|---|---|
| zero pixels, tall forest | 39.8 % | 17.1 % | 9.9 % | 13.9–14.0 % |
| cover, tall forest (median) | 0.60 | 0.86 | 0.95 | 0.90 |
| p90 / SD / mean, forest cells (m) | 25.9 / 10.2 / 9.3 | 28.9 / 9.4 / 15.9 | 30.4 / 7.9 / 20.1 | 20.6–21.1 / 6.7–6.9 / 11.3–12.3 |
| cover, open land (mean) | 0.085 | 0.118 | 0.145 | |

3 m removes the speckle and keeps clearings, forest lanes and open land; 5 m smears crowns and saturates cover. **3 m
is used** (`GAP_FILL_SIZE`, GeoTIFF tag `gap_fill_closing_m`). Cover and the 0.3 forest threshold become comparable
with England; p90 and mean height in forest rise by about 3–6 m (filled pixels take the height of the surrounding
crown), so absolute heights are not exactly comparable with England.

### 4.2 Tall structures: 60 m maximum

The EA VOM holds vegetation only. DSM − DTM also contains cranes, wind turbines, masts and chimneys (up to 191 m).
Heights above 60 m are set to 0 (`MAX_VEGETATION_HEIGHT`, tag `max_vegetation_height_m`); the tallest Belgian trees
are about 50 m. On three test blocks this affected 0.014 % of pixels or less, most of them on buildings already.
Known limitation: industrial structures of 50–60 m in the ports (Antwerp, Zeebrugge, Ghent) remain; a few hundred
30 m cells in Belgium have p90 > 50 m.

## 5. Results and checks

Checked with `scripts/check_metrics.py` (grid, tags, value ranges, seams, map) and `scripts/check_borders.py` (steps at
the regional borders).

| 30 m cells | Flanders | Brussels | Wallonia | Belgium |
|---|---|---|---|---|
| cells with data | 15.2 M | 0.18 M | 18.8 M | 34.1 M (30,718 km²) |
| forest cells (cover ≥ 0.3) | 20.2 % | 30.6 % | 33.7 % | 27.7 % |
| in forest cells, median p90 / SD / mean | 17.4 / 5.5 / 9.1 m | 16.8 / 6.3 / 7.0 m | 21.0 / 5.5 / 13.5 m | 19.6 / 5.5 / 11.8 m |

- Forest cells in Flanders (20 %) exceed the official forest area (about 11 %): 30 m cells with ≥ 30 % tree cover also
  include tree rows, parks and wooded gardens. Wallonia (34 %) is close to its forest share of about a third.
- **Seams:** the difference between neighbouring cells across processing tiles and download blocks equals that
  elsewhere (Belgium: ratios 1.00–1.02 at the 6 km tiles, 0.94–1.00 at the 10 km download blocks), so there are
  no edge effects.
- **Same results per region and for Belgium:** tiles inside one region are identical in single-region and all-Belgium
  runs.
- **Regional borders** (bands of 500 m on either side, forest cells): Flanders → Brussels (same survey) cover −0.02,
  p90 −0.4 m; Flanders → Wallonia cover −0.13, p90 −0.9 m, mean −2.3 m (see section 6).
- **Wallonia against DHMV II, Grez-Doiceau** (a 6 km box both surveys cover): correlation per cell 0.98–0.99 (cover,
  p90, mean) and 0.96–0.98 (SD); Wallonia is slightly higher (cover +0.04, p90 +1.2–1.4 m) because its DSM represents
  crowns slightly more fully (ground and roofs agree within 4 cm). This offset is accepted.

## 6. Data quality in Wallonia

### 6.1 Both surveys along the Flemish–Walloon border

DHMV II also covers a strip of Wallonia. Within 2 km of Flanders (603 km², `scripts/download_border_strip.py`), both
surveys were processed with the same rules and compared on 655,432 complete 30 m cells
(`scripts/compare_border_strip.py`):

| forest cells, Wallonia minus DHMV II | cover | p90 | mean height |
|---|---|---|---|
| whole strip | −0.04 | +0.7 m | +0.4 m |
| west, Hainaut (x 40–120 km) | −0.10 to −0.34 | −1.9 to −5.0 m | −1.3 to −7.0 m |
| centre and east (x 140–260 km) | +0.04 to +0.06 | +1.3 to +4.8 m | +1.3 to +4.6 m |

The difference is not one offset: it changes sign along the border. In western Hainaut the Walloon DSM lies on the
ground in deciduous forest (10 % of pixels above 2.5 m against 58 % in DHMV II; DTM equal within 4 cm): leaf-off
crowns were hardly recorded by those flights, and gap filling cannot repair this. Point density does not explain it;
flight dates per block are the likely cause (SPW publishes no flight blocks for this survey).

### 6.2 Extent inside Wallonia

Inside Wallonia there is no second survey, so the SPW forest mask 2018 (independent of both LiDAR surveys) stands in
for it (`scripts/download_forest_mask.py`, `scripts/check_wallonia_quality.py`). Measure: in cells entirely forest in
the mask, the share that the survey finds as forest (cover ≥ 0.3).

| in mask-forest cells | DHMV II | Wallonia 2013–14 |
|---|---|---|
| border strip, all sections | 95 % found | 90 % |
| border strip, section x 100–120 km (Hainaut) | 99 % | 64 % |
| per 6 km tile (502 tiles): median / 10 % / 25 % | – | 92 % / 36 % / 84 % |

The tiles fall into two groups: most find 80–100 % of mask forest (as good as DHMV II), about 90 tiles find 10–60 %,
with straight edges between them: flight blocks. Outside these blocks the survey is good, including broadleaf-rich
regions (Condroz, Famenne), so the problem is the flights, not broadleaf forest as such.

### 6.3 Quality layer

`scripts/make_quality_mask.py` draws the zone from the data as three rectangles (EPSG:31370), with edges on whole
kilometres where the share of mask forest found jumps, clipped to Wallonia:

| rectangle (km) | area |
|---|---|
| x < 162, y 76–138 | Hainaut south of Tournai and Ath (Mons, Charleroi, Thuin) |
| x 162–194, y 76–98 | Entre-Sambre-et-Meuse west of the Meuse |
| x 84–118, y ≥ 158 | the strip around Flobecq, Lessines and Enghien |

| code | meaning | cells (area) | in SPW mask forest, found as forest |
|---|---|---|---|
| 0 | no data | – | – |
| 1 | DHMV II (Flanders, Brussels) | 15.35 M (13,816 km²) | 95 % (border strip) |
| 2 | Wallonia 2013–14, reliable | 15.18 M (13,660 km²) | 93.1 % |
| 3 | Wallonia 2013–14, leaf-off broadleaf forest missed | 3.60 M (3,241 km²) | 44.2 % |

The zone covers 19 % of Wallonia and 16 % of the Walloon forest. A cell takes the region its centre lies in. Replacing
the zone with the 2021–2022 survey remains possible in a later version.

## 7. Notes for using the layers

- **Code 3 cells:** cover, mean height and p90 are far too low in broadleaf forest there. Leave them out when fitting a
  model (presence and background points), and mark or clip the zone in predictions. A forest mask from cover ≥ 0.3
  also misses forest there.
- **Heights are not exactly comparable with England:** after gap filling p90 and mean in forest are 3–6 m higher than
  without; Walloon heights are about 1 m above Flemish ones on the same ground. A Belgian model fitted on these layers
  is unaffected by the first point; a transfer of the England model is affected by both, and also needs England's PCA
  centring and scaling.
- **Reprojection:** cover and the 0.3 threshold are sensitive to bilinear resampling; producing the layers on the
  model's own grid (`GRID_ORIGIN`, CRS) avoids a second resampling.
- **Date:** the surveys date from 2012–2015; later forest change (felling, storm and bark-beetle damage) is not in the
  layers.

## Sources

[EA VOM](https://www.data.gov.uk/dataset/lidar-vegetation-object-model-vom) ·
[Digitaal Hoogtemodel Vlaanderen](https://www.vlaanderen.be/digitaal-vlaanderen/onze-diensten-en-platformen/earth-observation-data-science-eodas/het-digitaal-hoogtemodel) ·
[Géoportail de la Wallonie](https://geoportail.wallonie.be)
