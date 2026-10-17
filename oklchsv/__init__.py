"""OklchSV -- closed-form and spectral S/P ratio, with D65 adaptation.

This package contains only the S/P computation that supports
arXiv:2609.24819. The OklchSV mesopic transform itself is not included.
"""

from .sp import (YP, YPP, sp_closed_form, sp_spectral, fill_undefined,
                 load_cmf, load_scotopic_v)
from .adapt import (D65_XY, D65_XYZ, xy_to_XYZ, xyz_to_linear_srgb,
                    linear_srgb_to_xyz, bradford_matrix, adapt_to_d65)

__version__ = "0.1.0"

__all__ = [
    "YP", "YPP", "sp_closed_form", "sp_spectral", "fill_undefined",
    "load_cmf", "load_scotopic_v",
    "D65_XY", "D65_XYZ", "xy_to_XYZ", "xyz_to_linear_srgb",
    "linear_srgb_to_xyz", "bradford_matrix", "adapt_to_d65",
]
