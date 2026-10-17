# OklchSV — S/P ratio (v0.1)

[日本語](README.ja.md)

Code for the scotopic-to-photopic luminance ratio (S/P) used by OklchSV, a
mesopic display model. It accompanies the paper

> **When is a closed-form RGB->S/P ratio adequate? A hyperspectral characterization on natural scenes for mesopic display**
> arXiv:2609.24819 (CC BY 4.0) — <https://arxiv.org/abs/2609.24819>

This release contains the S/P part only:

- the **closed-form S/P** from linear sRGB (six published weights),
- the **spectral S/P** from radiance (V(λ) and V′(λ) integrals),
- the **D65 adaptation** the closed form needs (linear sRGB ⇄ XYZ, Bradford).

The OklchSV mesopic transform itself is **not** included.

## Install

```bash
pip install -e .              # numpy only
pip install -e ".[scripts]"   # + scipy, colour-science (for scripts/)
```

Python 3.10 or later.

## Use

An ordinary 8-bit sRGB image (Pillow is only used for reading the file here):

```python
import numpy as np
from PIL import Image
from oklchsv import srgb_to_linear, sp_closed_form

img = np.asarray(Image.open("photo.jpg").convert("RGB")) / 255.0   # encoded sRGB
rgb = srgb_to_linear(img)            # decode the gamma first
sp, ok = sp_closed_form(rgb)         # (H, W) S/P map; white and greys -> 1.0
```

$$
S/P = \frac{(R,G,B)\cdot Y'_p \,/\, \textstyle\sum Y'_p}{(R,G,B)\cdot Y_p},\qquad
Y_p = (0.2126,\ 0.7152,\ 0.0722),\quad Y'_p = (0.0184,\ 0.5673,\ 0.3498)
$$

The input must be **linear** sRGB (not gamma-encoded; use `srgb_to_linear`)
**adapted to D65** (paper condition C2). Ordinary photographs and museum images
that are white-balanced can be treated as D65-adapted. For a scene under another
illuminant, adapt first:

```python
from oklchsv import adapt_to_d65, linear_srgb_to_xyz, xyz_to_linear_srgb

rgb = xyz_to_linear_srgb(adapt_to_d65(linear_srgb_to_xyz(rgb), white_XYZ))
```

The spectral S/P (`sp_spectral`) needs the CIE tables, which are not bundled.
See [data/README.md](data/README.md).

## Scripts

| Script | Reproduces | Needs |
|---|---|---|
| `scripts/oklchsv_scotopic_reconcile.py` | Paper Sec. 3.2–3.3, Table 1: derivation of Y′ₚ from Gaussian-primary spectra, and the comparison with Maksimainen et al. | scipy, colour-science 0.4.7 |
| `scripts/hsi_sp_maps.py` | Closed-form vs spectral S/P maps on the Levada time-lapse HSI and the median difference (paper Stage 2) | scipy, CIE tables, Levada data |

`oklchsv_scotopic_reconcile.py` is a byte-identical copy of the arXiv ancillary file of the
same name (SHA-256 `5a0b01099f9ffc8b355b091f333a6a127054c84892f5d81f6fcc3e425da5c71a`).
The arXiv ancillary files are the reference for the paper; this copy is for convenience.

```bash
python scripts/oklchsv_scotopic_reconcile.py
python scripts/hsi_sp_maps.py                 # data in data/ (see data/README.md)
python scripts/hsi_sp_maps.py --save-maps output/
```

Expected median |closed-form − spectral| from `hsi_sp_maps.py`:

| Levada | 14:11 noon | 17:12 evening | 18:45 dusk |
|---|---|---|---|
| median \|ΔS/P\| | 0.0696 | 0.0657 | 0.0556 |

The paper's evaluation used Foster's `xyzbar.mat` for x̄ȳz̄ (see data/README.md).
The same four-digit values were also obtained with the CIE 1931 table as distributed
in colour-science 0.4.7. Other copies of the tables may differ in the last digit.

Where the spectral S/P is not defined (photopic luminance 0), the saved spectral map
takes the closed-form value. The medians above are over pixels where both are defined.

## What v0.1 does not reproduce

- Any rendering by the OklchSV transform (paper Fig. 5 is reproduced by the arXiv
  ancillary script with illustrative constants, not by this repository).
- The night renderings in the companion article (closed-form-driven vs spectral-driven
  appearance). They need the transform, which is not public.

## Tests

```bash
python -m unittest discover tests
```

Tests that need the CIE tables are skipped when the tables are not in `data/`.

## Citation

Please cite the paper (see [CITATION.cff](CITATION.cff)).

## License

MIT (code). The paper is CC BY 4.0. Data and CIE tables are not included; see
data/README.md for their sources and terms.
