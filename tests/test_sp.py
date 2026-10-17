"""Tests for oklchsv.sp. Run: python -m unittest discover tests"""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from oklchsv import sp  # noqa: E402

DATA = ROOT / "data"
WL = np.arange(400, 721, 10, dtype=float)
HAVE_TABLES = ((DATA / "cie_1931_2deg_xyz.csv").is_file()
               or (DATA / "Tutorial_HSI2RGB" / "xyzbar.mat").is_file()) \
    and (DATA / "cie_1951_scotopic_v.csv").is_file()


class TestClosedForm(unittest.TestCase):
    def test_weights(self):
        np.testing.assert_allclose(sp.YP, [0.2126, 0.7152, 0.0722])
        np.testing.assert_allclose(sp.YPP, [0.0184, 0.5673, 0.3498])
        self.assertAlmostEqual(sp.YP.sum(), 1.0, places=12)
        self.assertAlmostEqual(sp.YPP_SUM, 0.9355, places=12)

    def test_white_is_one(self):
        for k in (1.0, 0.18, 250.0):
            s, ok = sp.sp_closed_form([k, k, k])
            self.assertTrue(ok)
            self.assertAlmostEqual(float(s), 1.0, places=12)

    def test_primaries(self):
        s, _ = sp.sp_closed_form(np.eye(3))
        expect = (sp.YPP / 0.9355) / sp.YP
        np.testing.assert_allclose(s, expect, rtol=1e-12)
        # blue > green > red in S/P
        self.assertTrue(s[2] > s[1] > s[0])

    def test_scale_invariance(self):
        rgb = np.array([0.3, 0.5, 0.1])
        s1, _ = sp.sp_closed_form(rgb)
        s2, _ = sp.sp_closed_form(7.0 * rgb)
        self.assertAlmostEqual(float(s1), float(s2), places=12)

    def test_undefined(self):
        s, ok = sp.sp_closed_form([[0, 0, 0], [1, 1, 1]])
        self.assertFalse(ok[0])
        self.assertTrue(np.isnan(s[0]))
        self.assertTrue(ok[1])


class TestSpectral(unittest.TestCase):
    """Definition checks with synthetic V, V' (no CIE tables needed)."""

    def setUp(self):
        rng = np.random.default_rng(0)
        self.V = rng.uniform(0.1, 1.0, 5)
        self.Vp = rng.uniform(0.1, 1.0, 5)
        self.white = rng.uniform(0.5, 1.5, 5)

    def test_white_is_one(self):
        s = sp.sp_spectral(3.0 * self.white, self.V, self.Vp, self.white)
        self.assertAlmostEqual(float(s), 1.0, places=12)

    def test_definition(self):
        L = np.array([1.0, 2.0, 0.5, 0.0, 3.0])
        s = sp.sp_spectral(L, self.V, self.Vp, self.white)
        expect = ((L @ self.Vp) / (self.white @ self.Vp)) / \
                 ((L @ self.V) / (self.white @ self.V))
        self.assertAlmostEqual(float(s), float(expect), places=12)

    def test_zero_luminance_not_finite(self):
        s = sp.sp_spectral(np.zeros((2, 5)), self.V, self.Vp, self.white)
        self.assertFalse(np.isfinite(s).any())

    def test_fill_undefined(self):
        prim = np.array([1.2, np.nan, np.inf, np.nan])
        fb = np.array([9.0, 0.8, 0.9, np.nan])
        out, n = sp.fill_undefined(prim, fb)
        np.testing.assert_allclose(out, [1.2, 0.8, 0.9, 1.0])
        self.assertEqual(n, 3)


class TestTables(unittest.TestCase):
    def test_missing_table_message(self):
        with self.assertRaises(FileNotFoundError) as cm:
            sp.load_scotopic_v(ROOT / "no_such_table.csv", WL)
        self.assertIn("data/README.md", str(cm.exception))

    @unittest.skipUnless(HAVE_TABLES, "CIE tables not in data/ (see data/README.md)")
    def test_tables_on_levada_grid(self):
        cmf_path = DATA / "Tutorial_HSI2RGB" / "xyzbar.mat"
        if not cmf_path.is_file():
            cmf_path = DATA / "cie_1931_2deg_xyz.csv"
        cmf = sp.load_cmf(cmf_path, WL)
        Vp = sp.load_scotopic_v(DATA / "cie_1951_scotopic_v.csv", WL)
        self.assertEqual(cmf.shape, (33, 3))
        self.assertEqual(Vp.shape, (33,))
        self.assertEqual(WL[np.argmax(cmf[:, 1])], 560.0)   # V peaks at 555 nm
        self.assertEqual(WL[np.argmax(Vp)], 510.0)          # V' peaks at 507 nm
        # equal-energy spectrum vs itself as white -> 1
        E = np.ones(33)
        self.assertAlmostEqual(float(sp.sp_spectral(E, cmf[:, 1], Vp, E)), 1.0)


if __name__ == "__main__":
    unittest.main()
