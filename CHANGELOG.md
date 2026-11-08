# Changelog

## Unreleased
* Documentation: the difference in input with England (DSM − DTM approximating the EA VOM, and why the VOM
  approach was not applied) is stated in the README, `docs/methods.md` (section 4.0) and the report (section 5.1);
  download link and data licence (CC BY 4.0) in the README; report section 1.4 is advice to the modeller.
  No change to the code or the layers. The v0.1.0 release files were made with the tagged code and are identical
  in every cell to the reference run.

## 0.1.0 (2026-11-08)
First version: the four 30 m forest-structure layers for all of Belgium and the quality layer.

* Canopy height from DSM − DTM encoded like the EA Vegetation Object Model (< 2.5 m and buildings → 0), with two
  deviations from England: 3 m gap filling of leaf-off crowns and a 60 m maximum height (`docs/methods.md`).
* Inputs: DHMV II (Flanders, Brussels; 2013–2015) and the Walloon survey 2013–2014.
* Quality layer: source survey per cell and the Walloon zone (code 3) where leaf-off broadleaf forest is missed.
* Grid EPSG:31370, cells on multiples of 30 m (provisional).
* Results identical to the reference run of 2026-11-08 (`metrics_30m_belgium.tif`); the script name recorded in the
  output tags is now `belgium-lidar-forest-metrics/scripts/...`.
