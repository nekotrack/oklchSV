# Changelog

## 0.1.0

First public release.

- `oklchsv.sp`: closed-form S/P (six published weights) and spectral S/P
  (V and V′ integrals, white-normalised); loaders for the CIE tables
  (not bundled).
- `oklchsv.adapt`: linear sRGB ⇄ XYZ and Bradford adaptation to D65.
- `scripts/oklchsv_scotopic_reconcile.py`: copy of the arXiv ancillary file
  (paper Sec. 3.2).
- `scripts/hsi_sp_maps.py`: closed-form vs spectral S/P maps on the Levada
  time-lapse HSI.
