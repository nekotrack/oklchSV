# data/

No data are included in this repository. Put the files below here yourself.

## 1. CIE tables (needed for the spectral S/P)

The colour-matching functions and the luminous-efficiency functions are not bundled.
Download them from the CIE or CVRL and save them under these names:

| File (in `data/`) | Content | Columns | Source |
|---|---|---|---|
| `cie_1931_2deg_xyz.csv` | CIE 1931 2° colour-matching functions | `wavelength_nm, xbar, ybar, zbar` | CIE data tables (<https://cie.co.at/data-tables>) or CVRL (<http://www.cvrl.org/>) |
| `cie_1951_scotopic_v.csv` | CIE 1951 scotopic luminous efficiency V′(λ) | `wavelength_nm, V'` | same |

Format:

- Plain CSV, one wavelength per row, wavelength in nm in the first column.
  Header and comment lines are skipped; empty cells are read as 0.
- Any spacing works (1 nm or 5 nm). The tables are picked at, or linearly
  interpolated to, the grid that is needed (400:10:720 nm for Levada). They must
  cover that range.
- The photopic V(λ) is taken from ȳ (column `ybar`), as in the paper.

If a table is missing, the code stops and tells you which file it looked for.
The tests that need the tables are skipped.

### Exact reproduction of the paper's numbers

The paper's evaluation took x̄ȳz̄ from `xyzbar.mat` in Foster's HSI→RGB tutorial
(400:10:720 nm). If you place it at `data/Tutorial_HSI2RGB/xyzbar.mat`,
`scripts/hsi_sp_maps.py` uses it in preference to `cie_1931_2deg_xyz.csv`
(reading it needs SciPy). The tutorial package is available from:

- <https://personalpages.manchester.ac.uk/staff/d.h.foster/Hyperspectral_images_of_natural_scenes_02.html>

## 2. Levada time-lapse HSI (for `scripts/hsi_sp_maps.py`)

Foster, Amano & Nascimento (2016), *Vision Research* 120, 45–60.
Time-lapse hyperspectral radiance images, 2015 set:

- <https://personalpages.manchester.ac.uk/staff/d.h.foster/Time-Lapse_HSIs/Time-Lapse_HSIs_2015.html>

These data are for personal use only and must not be redistributed. Extract the
three timestamps so that each folder contains its `.mat` file:

```
data/
├─ cie_1931_2deg_xyz.csv
├─ cie_1951_scotopic_v.csv
├─ Tutorial_HSI2RGB/xyzbar.mat     (optional, for exact reproduction)
└─ levada/
   ├─ levada_1411/levada_1411.mat
   ├─ levada_1712/levada_1712.mat
   └─ levada_1845/levada_1845.mat
```

Each cube is 1024 × 1344 × 33 (400:10:720 nm), spectral radiance in
W m⁻² sr⁻¹ nm⁻¹. One cube takes about 360 MB in memory as float64.
