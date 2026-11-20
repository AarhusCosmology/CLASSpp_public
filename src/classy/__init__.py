"""classy, the Python interface to CLASS++.

`Class` computes a cosmology from a dict of CLASS input parameters and answers questions about
it. The compiled core (`_core`) holds the C++ Cosmology and offers one method per C++ call;
everything else in this package is plain Python.
"""
from ._errors import CosmoError, CosmoSevereError, CosmoComputationError
from ._core import __version__, PyCosmology
from ._bessel import (hyperspherical_bessel_direct, hyperspherical_bessel_interpolate,
                    hyperspherical_bessel_wkb)
from ._class import Class

__all__ = ['Class', 'PyCosmology', 'CosmoError', 'CosmoSevereError', 'CosmoComputationError', '__version__',
           'hyperspherical_bessel_direct', 'hyperspherical_bessel_interpolate',
           'hyperspherical_bessel_wkb']
