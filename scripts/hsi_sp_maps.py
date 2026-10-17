#!/usr/bin/env python3
"""Closed-form vs spectral S/P maps on Foster's Levada time-lapse HSI.

For each timestamp:
  1. radiance (H, W, 33; 400:10:720 nm) -> photopic luminance [cd/m^2]
  2. scene white = mean radiance of the pixels with luminance in [p97, p99]
     (highlight band)
  3. spectral S/P: V(lambda) = ybar, V'(lambda) = CIE 1951, normalised by the
     scene white
  4. closed-form S/P: radiance -> XYZ -> Bradford (scene white -> D65)
     -> linear sRGB -> six-weight formula
  5. median |closed-form - spectral| over pixels where both are defined
     (the statistic reported for Stage 2 of arXiv:2609.24819)

Where the spectral S/P is not defined (photopic luminance 0), the saved
spectral map takes the closed-form value (and 1.0 where neither is defined);
the count is printed as n_filled. The median above does not use filled pixels.

The data are not included (personal use only); see data/README.md.

Usage:
  python scripts/hsi_sp_maps.py
  python scripts/hsi_sp_maps.py --data-dir path/to/data --save-maps output/
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from oklchsv import (adapt_to_d65, fill_undefined, load_cmf,  # noqa: E402
                     load_scotopic_v, sp_closed_form, sp_spectral,
                     xyz_to_linear_srgb)

WL = np.arange(400, 721, 10, dtype=float)     # Levada grid, 33 bands
KM = 683.0                                    # lm/W
DLAM = 10.0                                   # nm
TIMES = [1411, 1712, 1845]
TIME_LABEL = {1411: "14:11 noon", 1712: "17:12 evening", 1845: "18:45 dusk"}


def find_cmf(data_dir):
    """Prefer Foster's xyzbar.mat (the paper's table), else the CIE CSV."""
    for p in (data_dir / "Tutorial_HSI2RGB" / "xyzbar.mat",
              data_dir / "cie_1931_2deg_xyz.csv"):
        if p.is_file():
            return p
    return data_dir / "cie_1931_2deg_xyz.csv"     # load_cmf raises a clear error


def load_cube(data_dir, t):
    """First 3-D (., ., 33) variable of the .mat under levada/levada_<t>/."""
    from scipy.io import loadmat
    tdir = data_dir / "levada" / f"levada_{t}"
    mats = sorted(tdir.rglob("*.mat")) if tdir.is_dir() else []
    if not mats:
        raise FileNotFoundError(
            f"no .mat under {tdir}. See data/README.md for the Levada files.")
    m = loadmat(mats[0])
    for k, v in m.items():
        if not k.startswith("__") and np.ndim(v) == 3 and v.shape[-1] == WL.size:
            return np.asarray(v, dtype=float)
    raise ValueError(f"{mats[0]}: no (H, W, {WL.size}) variable")


def highlight_white(radiance_flat, L_cd):
    """Mean spectrum of the pixels with luminance in [p97, p99]."""
    fin = np.isfinite(L_cd)
    lo, hi = np.percentile(L_cd[fin], [97, 99])
    sel = fin & (L_cd >= lo) & (L_cd <= hi)
    return radiance_flat[sel].mean(axis=0), int(sel.sum())


def process(data_dir, t, cmf, Vp):
    rad = load_cube(data_dir, t)
    H, W, _ = rad.shape
    flat = rad.reshape(-1, WL.size)
    del rad
    V = cmf[:, 1]

    L_cd = KM * (flat @ V) * DLAM
    white_spd, n_white = highlight_white(flat, L_cd)
    XYZ_white = white_spd @ cmf
    XYZ_white = XYZ_white / XYZ_white[1]

    rgb = xyz_to_linear_srgb(adapt_to_d65(flat @ cmf, XYZ_white))
    sp_cf, ok_cf = sp_closed_form(rgb)
    sp_sp = sp_spectral(flat, V, Vp, white_spd)

    valid = ok_cf & np.isfinite(sp_cf) & np.isfinite(sp_sp)
    diff = (sp_cf - sp_sp)[valid]
    sp_sp_filled, n_filled = fill_undefined(sp_sp, sp_cf)
    return {
        "t": t,
        "median_abs_diff": float(np.median(np.abs(diff))),
        "median_signed_diff": float(np.median(diff)),
        "n_valid": int(valid.sum()),
        "n_pixels": int(valid.size),
        "n_filled": n_filled,
        "n_white": n_white,
        "maps": (sp_cf.reshape(H, W), sp_sp_filled.reshape(H, W)),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data")
    ap.add_argument("--times", type=int, nargs="+", default=TIMES)
    ap.add_argument("--save-maps", type=Path, default=None,
                    help="directory for levada_<t>_sp_maps.npz (closed, spectral)")
    args = ap.parse_args(argv)

    data_dir = args.data_dir.resolve()
    cmf_path = find_cmf(data_dir)
    cmf = load_cmf(cmf_path, WL)
    Vp = load_scotopic_v(data_dir / "cie_1951_scotopic_v.csv", WL)
    print(f"CMF: {cmf_path.name}")

    print(f"{'time':<15}{'median |diff|':>14}{'n_valid':>10}{'n_filled':>10}")
    for t in args.times:
        r = process(data_dir, t, cmf, Vp)
        print(f"{TIME_LABEL.get(t, str(t)):<15}{r['median_abs_diff']:>14.4f}"
              f"{r['n_valid']:>10d}{r['n_filled']:>10d}")
        if args.save_maps is not None:
            args.save_maps.mkdir(parents=True, exist_ok=True)
            out = args.save_maps / f"levada_{t}_sp_maps.npz"
            np.savez_compressed(out, closed_form=r["maps"][0],
                                spectral=r["maps"][1])
            print(f"  saved {out}")


if __name__ == "__main__":
    main()
