"""Scotopic-to-photopic luminance ratio (S/P).

Two definitions, both white-normalised so that S/P = 1 at the reference white:

* ``sp_closed_form`` -- the closed-form S/P from linear sRGB, using the six
  published scalar weights (arXiv:2609.24819, Sec. 3.1).
* ``sp_spectral`` -- the spectral S/P from a radiance spectrum, integrating
  V(lambda) and V'(lambda) and normalising by a white spectrum.

The closed form is meant to be evaluated on *adapted* linear sRGB (D65 white);
see :mod:`oklchsv.adapt` (paper condition C2).

The colour-matching and luminous-efficiency tables are not bundled. Use
``load_cmf`` / ``load_scotopic_v`` with files obtained as described in
``data/README.md``.
"""

from pathlib import Path

import numpy as np

# Photopic weights: the Y row of the sRGB (D65) RGB -> XYZ matrix. Sum = 1.
YP = np.array([0.2126, 0.7152, 0.0722])
# Scotopic weights: CIE 1951 V'(lambda) projected onto Gaussian-primary
# sRGB spectra (paper Sec. 3.2). Sum = 0.9355; divided out below.
YPP = np.array([0.0184, 0.5673, 0.3498])
YPP_SUM = float(YPP.sum())

DATA_README = "data/README.md"


def sp_closed_form(rgb, eps=1e-9):
    """Closed-form S/P of linear sRGB.

    S/P = ((rgb . YPP) / sum(YPP)) / (rgb . YP), which is 1 at (1, 1, 1).

    Parameters
    ----------
    rgb : array_like, shape (..., 3)
        Linear (not gamma-encoded) sRGB, adapted to D65.
    eps : float
        Pixels whose photopic denominator is <= eps are not defined.

    Returns
    -------
    sp : ndarray, shape (...)
        S/P, NaN where undefined.
    ok : ndarray of bool, shape (...)
        True where S/P is defined.
    """
    rgb = np.asarray(rgb, dtype=float)
    num = (rgb @ YPP) / YPP_SUM
    den = rgb @ YP
    ok = den > eps
    sp = np.full(den.shape, np.nan)
    sp[ok] = num[ok] / den[ok]
    return sp, ok


def sp_spectral(radiance, V, Vp, white_spd):
    """Spectral S/P, white-normalised.

    S/P = (Ys / Ys_white) / (Yp / Yp_white), with Yp = L . V and
    Ys = L . V'. The band width and the luminous efficacies cancel.

    Parameters
    ----------
    radiance : array_like, shape (..., n)
        Spectral radiance (any linear unit) on the same wavelength grid as V.
    V, Vp : array_like, shape (n,)
        Photopic V(lambda) and scotopic V'(lambda) on that grid.
    white_spd : array_like, shape (n,)
        Spectrum of the reference white.

    Returns
    -------
    sp : ndarray, shape (...)
        S/P. Non-finite where the photopic luminance is 0 (division by zero
        is not raised; check with ``np.isfinite``).
    """
    radiance = np.asarray(radiance, dtype=float)
    V = np.asarray(V, dtype=float)
    Vp = np.asarray(Vp, dtype=float)
    white_spd = np.asarray(white_spd, dtype=float)
    Yp_w = float(white_spd @ V)
    Ys_w = float(white_spd @ Vp)
    with np.errstate(divide="ignore", invalid="ignore"):
        return ((radiance @ Vp) / Ys_w) / ((radiance @ V) / Yp_w)


def fill_undefined(sp_primary, sp_fallback, default=1.0):
    """Replace non-finite values of ``sp_primary`` by ``sp_fallback``.

    Used to make a complete spectral S/P map: where the spectral S/P is not
    defined, take the closed-form value; where that is also undefined, take
    ``default`` (S/P of the white).

    Returns
    -------
    filled : ndarray
    n_filled : int
        Number of values taken from ``sp_fallback`` or ``default``.
    """
    sp_primary = np.asarray(sp_primary, dtype=float)
    bad = ~np.isfinite(sp_primary)
    filled = np.where(bad, sp_fallback, sp_primary)
    filled = np.where(np.isfinite(filled), filled, default)
    return filled, int(bad.sum())


# ---------------------------------------------------------------------------
#  table loading (files are not bundled; see data/README.md)
# ---------------------------------------------------------------------------
def _require(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"table not found: {path}\n"
            f"The CIE tables are not bundled with this repository. "
            f"See {DATA_README} for where to get them and how to name them.")
    return path


def _read_csv_table(path, ncols):
    """Read 'wavelength, v1[, v2, ...]' rows; header lines are skipped."""
    rows = []
    with open(_require(path), encoding="utf-8-sig") as f:
        for line in f:
            parts = [p.strip() for p in line.replace(";", ",").split(",")]
            if len(parts) < ncols + 1:
                continue
            try:
                vals = [float(p) if p else 0.0 for p in parts[:ncols + 1]]
            except ValueError:
                continue                    # header or comment line
            rows.append(vals)
    if not rows:
        raise ValueError(f"no numeric rows in {path}")
    a = np.array(rows)
    return a[:, 0], a[:, 1:]


def _on_grid(wl_src, values, wavelengths):
    """Pick (or linearly interpolate) the table onto ``wavelengths``."""
    wavelengths = np.asarray(wavelengths, dtype=float)
    if wavelengths.min() < wl_src.min() or wavelengths.max() > wl_src.max():
        raise ValueError(
            f"table covers {wl_src.min():g}-{wl_src.max():g} nm, "
            f"requested {wavelengths.min():g}-{wavelengths.max():g} nm")
    return np.stack([np.interp(wavelengths, wl_src, values[:, j])
                     for j in range(values.shape[1])], axis=1)


def load_cmf(path, wavelengths):
    """CIE 1931 2-degree colour-matching functions on ``wavelengths``.

    ``path`` is either a CSV ('wavelength, xbar, ybar, zbar') or the
    ``xyzbar.mat`` of Foster's HSI tutorial (33 x 3, 400:10:720 nm; reading
    it needs SciPy).

    Returns
    -------
    cmf : ndarray, shape (n, 3)
        Columns xbar, ybar (= photopic V), zbar.
    """
    path = _require(path)
    wavelengths = np.asarray(wavelengths, dtype=float)
    if path.suffix.lower() == ".mat":
        from scipy.io import loadmat
        cmf = np.asarray(loadmat(path)["xyzbar"], dtype=float)
        grid = np.arange(400, 721, 10, dtype=float)
        if cmf.shape != (grid.size, 3):
            raise ValueError(f"unexpected xyzbar shape {cmf.shape}")
        return _on_grid(grid, cmf, wavelengths)
    wl, vals = _read_csv_table(path, 3)
    return _on_grid(wl, vals, wavelengths)


def load_scotopic_v(path, wavelengths):
    """CIE 1951 scotopic luminous efficiency V'(lambda) on ``wavelengths``.

    ``path`` is a CSV ('wavelength, V'').
    """
    wl, vals = _read_csv_table(path, 1)
    return _on_grid(wl, vals, wavelengths)[:, 0]
