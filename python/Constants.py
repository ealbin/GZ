#!/usr/bin/env python

"""Wrapper for physical constants.
"""

import sys

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"


# physical constants and conversions:
__c__           = 299792458.    # speed of light [meters / second]
__ev_per_amu__  = 931.5e6       # [eV/c**2] mass of 1 atomic mass unit [amu] 
__g_per_t__     = 1e4           # magnetic field [Gauss] in 1 [Tesla]
__m_per_AU__    = 149597870700. # [meters / astronomical unit]
__r_sun_m__     = 695508000.    # radius of the Sun [meters]
__r_earth_m__   = 6378100.      # radius of Earth [meters]


def lightspeed_m_per_s():
    """ Returns the speed of light in [meters / second]
    """
    return __c__

def amu2eV(amu):
    """ Converts atomic mass units [amu] into [eV/c**2] mass
    """
    return amu * __ev_per_amu__

def Gauss2Tesla(gauss):
    """ Converts magnetic field [Gauss] into [Tesla]
    """
    return gauss / __g_per_t__

def Tesla2Gauss(tesla):
    """ Converts magnetic field [Tesla] into [Gauss]
    """
    return tesla * __g_per_t__

def meters2AU(meters):
    """ Converts a distance "meters" into its equivalent in [Astronomical Units]
    """
    return meters / __m_per_AU__

def AU2meters(AU):
    """ Converts a distance in "AU" into its equivalent in [meters]
    """
    return AU * __m_per_AU__

def RadiusMeters(who):
    """ Returns the radius of either the "earth" or "sun" in [meters]
    """
    if who.lower() == 'sun':
        return __r_sun_m__
    elif who.lower() == 'earth':
        return __r_earth_m__
    else:
        print 'earth or sun'
        sys.exit(1)
    
def RadiusAU(who):
    """ Returns the radius of either the "earth" or "sun" in [Astronomical Units]
    """
    if who == 'sun':
        return meters2AU(__r_sun_m__)
    elif who == 'earth':
        return meters2AU(__r_earth_m__)
    else:
        print 'earth or sun'
        sys.exit(1)
