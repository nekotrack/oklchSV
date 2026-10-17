"""Linear sRGB <-> XYZ and Bradford chromatic adaptation to D65.

The closed-form S/P assumes linear sRGB adapted to the D65 white
(arXiv:2609.24819, condition C2). For a scene under another illuminant,
adapt its XYZ to D65 first (``adapt_to_d65``), then convert to linear sRGB.

NumPy only. The Bradford step is a von Kries scaling in the Bradford cone
space, applied to XYZ as given (whites are not renormalised).
"""

import numpy as np

# XYZ -> linear sRGB, 4-digit coefficients of IEC 61966-2-1 (as in the
# paper's evaluation).
M_XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986],
                          [-0.9689, 1.8758, 0.0414],
                          [0.0557, -0.2040, 1.0570]])
M_SRGB_TO_XYZ = np.linalg.inv(M_XYZ_TO_SRGB)

M_BRADFORD = np.array([[0.8951, 0.2664, -0.1614],
                       [-0.7502, 1.7135, 0.0367],
                       [0.0389, -0.0685, 1.0296]])

# CIE 1931 2-degree chromaticity of D65.
D65_XY = (0.3127, 0.3290)


def xy_to_XYZ(xy, Y=1.0):
    """Chromaticity (x, y) -> XYZ with luminance ``Y``."""
    x, y = float(xy[0]), float(xy[1])
    return np.array([x / y * Y, Y, (1.0 - x - y) / y * Y])


D65_XYZ = xy_to_XYZ(D65_XY)


def xyz_to_linear_srgb(XYZ):
    """XYZ (..., 3) -> linear sRGB (..., 3). No clipping."""
    return np.asarray(XYZ, dtype=float) @ M_XYZ_TO_SRGB.T


def linear_srgb_to_xyz(rgb):
    """Linear sRGB (..., 3) -> XYZ (..., 3)."""
    return np.asarray(rgb, dtype=float) @ M_SRGB_TO_XYZ.T


def bradford_matrix(white_src, white_dst=D65_XYZ):
    """3x3 matrix adapting XYZ from ``white_src`` to ``white_dst`` (XYZ)."""
    lms_src = M_BRADFORD @ np.asarray(white_src, dtype=float)
    lms_dst = M_BRADFORD @ np.asarray(white_dst, dtype=float)
    return np.linalg.inv(M_BRADFORD) @ np.diag(lms_dst / lms_src) @ M_BRADFORD


def adapt_to_d65(XYZ, white_XYZ):
    """Bradford-adapt XYZ (..., 3) from the scene white ``white_XYZ`` to D65."""
    M = bradford_matrix(white_XYZ, D65_XYZ)
    return np.asarray(XYZ, dtype=float) @ M.T
