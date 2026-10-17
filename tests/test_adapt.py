"""Tests for oklchsv.adapt. Run: python -m unittest discover tests"""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from oklchsv import adapt  # noqa: E402

try:
    import colour
except ImportError:          # optional cross-check
    colour = None

ILL_A_XY = (0.44757, 0.40745)


class TestAdapt(unittest.TestCase):
    def test_d65_white_does_not_move(self):
        M = adapt.bradford_matrix(adapt.D65_XYZ, adapt.D65_XYZ)
        np.testing.assert_allclose(M, np.eye(3), atol=1e-12)
        xyz = np.array([[0.2, 0.3, 0.1], [0.9, 1.0, 1.1]])
        np.testing.assert_allclose(adapt.adapt_to_d65(xyz, adapt.D65_XYZ), xyz,
                                   atol=1e-12)

    def test_source_white_maps_to_d65(self):
        w = adapt.xy_to_XYZ(ILL_A_XY)
        np.testing.assert_allclose(adapt.adapt_to_d65(w, w), adapt.D65_XYZ,
                                   atol=1e-12)

    def test_srgb_round_trip(self):
        rgb = np.random.default_rng(1).uniform(0, 1, (10, 3))
        back = adapt.xyz_to_linear_srgb(adapt.linear_srgb_to_xyz(rgb))
        np.testing.assert_allclose(back, rgb, atol=1e-12)

    def test_srgb_white_is_d65(self):
        # 4-digit matrix: agreement to ~1e-4
        np.testing.assert_allclose(adapt.linear_srgb_to_xyz([1, 1, 1]),
                                   adapt.D65_XYZ, atol=2e-4)

    def test_d65_white_is_srgb_white_after_adaptation(self):
        w = adapt.xy_to_XYZ(ILL_A_XY)
        rgb = adapt.xyz_to_linear_srgb(adapt.adapt_to_d65(w, w))
        np.testing.assert_allclose(rgb, [1, 1, 1], atol=2e-4)

    @unittest.skipIf(colour is None, "colour-science not installed")
    def test_matches_colour_science(self):
        d65 = colour.xy_to_XYZ(
            colour.CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["D65"])
        np.testing.assert_allclose(adapt.D65_XYZ, d65, atol=1e-12)
        w = adapt.xy_to_XYZ(ILL_A_XY)
        xyz = np.random.default_rng(2).uniform(0, 1, (20, 3))
        ref = colour.chromatic_adaptation(xyz, w, d65, method="Von Kries",
                                          transform="Bradford")
        np.testing.assert_allclose(adapt.adapt_to_d65(xyz, w), ref, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
