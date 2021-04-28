"""
Source: http://www.wave-scattering.com/drudefit.html
"""

from .universal_constants import *
from .vacuum_constants import *


def eV_to_Hz(f):
    f *= ELECTRON_CHARGE #Joules
    f /= PLANCK_CONSTANT #Hertz
    return f

def metal_dict(plasma_freq, damping):
    fp = eV_to_Hz(plasma_freq)
    g = eV_to_Hz(damping)
    sigma = 2*PI*(fp**2)/g * VACUUM_PERMITTIVITY
    return {"conductivity":sigma, "plasma_frequency":fp, "damping":g}

GOLD = metal_dict(8.55, 0.0184)
SILVER = metal_dict(9.6, 0.0228)
ALUMINUM = metal_dict(15.3, 0.5984)
SODIUM = metal_dict(5.71, 0.0276)
POTASIUM = (3.72, 0.0184)
