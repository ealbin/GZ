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
__c__            = 299792458.    # speed of light [meters / second]
__ev_per_amu__   = 931.5e6       # [eV/c**2] mass of 1 atomic mass unit [amu] 
__g_per_t__      = 1e4           # magnetic field [Gauss] in 1 [Tesla]
__m_per_AU__     = 149597870700. # [meters / astronomical unit]
__r_sun_m__      = 695508000.    # radius of the Sun [meters]
__r_earth_m__    = 6378100.      # radius of Earth [meters]


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

def mass_number(proton_number):
    """ Returns an integer-rounded typical mass_number for a given proton_number.
    """
    if proton_number == 1:
        return 1 # Hydrogen
    elif proton_number == 2:
        return 4 # Helium
    elif proton_number == 3:
        return 7 # Lithium
    elif proton_number == 4:
        return 9 # Beryllium
    elif proton_number == 5:
        return 11 # Boron
    elif proton_number == 6:
        return 12 # Carbon
    elif proton_number == 7:
        return 14 # Nitrogen
    elif proton_number == 8:
        return 16 # Oxygen
    elif proton_number == 9:
        return 19 # Fluorine
    elif proton_number == 10:
        return 20 # Neon
    elif proton_number == 11:
        return 23 # Sodium
    elif proton_number == 12:
        return 24 # Magnesium
    elif proton_number == 13:
        return 27 # Alumininum
    elif proton_number == 14:
        return 28 # Silicon
    elif proton_number == 15:
        return 31 # Phosphorus
    elif proton_number == 16:
        return 32 # Sulfur
    elif proton_number == 17:
        return 35 # Chlorine
    elif proton_number == 18:
        return 40 # Argon
    elif proton_number == 19:
        return 39 # Potassium
    elif proton_number == 20:
        return 40 # Calcium
    elif proton_number == 21:
        return 45 # Scandium
    elif proton_number == 22:
        return 48 # Titanium
    elif proton_number == 23:
        return 51 # Vanadium
    elif proton_number == 24:
        return 52 # Chromium
    elif proton_number == 25:
        return 55 # Manganese
    elif proton_number == 26:
        return 56 # Iron
    elif proton_number == 27:
        return 59 # Cobalt
    elif proton_number == 28:
        return 59 # Nickel
    elif proton_number == 29:
        return 64 # Copper 
    elif proton_number == 30:
        return 65 # Zinc
    elif proton_number == 31:
        return 70 # Gallium
    elif proton_number == 32:
        return 73 # Germanium
    elif proton_number == 33:
        return 75 # Arsenic
    elif proton_number == 34:
        return 79 # Selenium
    elif proton_number == 35:
        return 80 # Bromine
    elif proton_number == 36:
        return 84 # Krypton
    elif proton_number == 37:
        return 85 # Rubidium
    elif proton_number == 38:
        return 88 # Strontium
    elif proton_number == 39:
        return 89 # Yttrium
    elif proton_number == 40:
        return 91 # Zirconium
    elif proton_number == 41:
        return 93 # Niobium
    elif proton_number == 42:
        return 96 # Molbdenum
    elif proton_number == 43:
        return 98 # Technium
    elif proton_number == 44:
        return 101 # Ruthenium
    elif proton_number == 45:
        return 103 # Rhodinium
    elif proton_number == 46:
        return 106 # Palladium
    elif proton_number == 47:
        return 108 # Silver
    elif proton_number == 48:
        return 112 # Cadmium
    elif proton_number == 49:
        return 115 # Indium
    elif proton_number == 50:
        return 119 # Tin
    elif proton_number == 51:
        return 122 # Antimony
    elif proton_number == 52:
        return 128 # Tellurium
    elif proton_number == 53:
        return 127 # Iodine
    elif proton_number == 54:
        return 131 # Xenon
    elif proton_number == 55:
        return 133 # Caesium
    elif proton_number == 56:
        return 137 # Barium
    elif proton_number == 57:
        return 139 # Lanthanum
    elif proton_number == 58:
        return 140 # Cerium
    elif proton_number == 59:
        return 141 # Praseodymium
    elif proton_number == 60:
        return 144 # Neodymium
    elif proton_number == 61:
        return 145 # Promethium
    elif proton_number == 62:
        return 150 # Samarium
    elif proton_number == 63:
        return 152 # Europium
    elif proton_number == 64:
        return 157 # Gadolinium
    elif proton_number == 65:
        return 159 # Terbium
    elif proton_number == 66:
        return 163 # Dysprosium
    elif proton_number == 67:
        return 165 # Holmium
    elif proton_number == 68:
        return 167 # Erbium
    elif proton_number == 69:
        return 169 # Thulium
    elif proton_number == 70:
        return 173 # Ytterbium
    elif proton_number == 71:
        return 175 # Lutetium
    elif proton_number == 72:
        return 178 # Hafnium
    elif proton_number == 73:
        return 181 # Tantalum
    elif proton_number == 74:
        return 184 # Tungsten
    elif proton_number == 75:
        return 186 # Rhenium
    elif proton_number == 76:
        return 190 # Osmium
    elif proton_number == 77:
        return 192 # Iridium
    elif proton_number == 78:
        return 195 # Platinum
    elif proton_number == 79:
        return 197 # Gold
    elif proton_number == 80:
        return 201 # Mercury
    elif proton_number == 81:
        return 204 # Thallium
    elif proton_number == 82:
        return 207 # Lead
    elif proton_number == 83:
        return 209 # Bismuth
    elif proton_number == 84:
        return 209 # Polonium
    elif proton_number == 85:
        return 210 # Astatine
    elif proton_number == 86:
        return 222 # Radon
    elif proton_number == 87:
        return 223 # Francium
    elif proton_number == 88:
        return 226 # Radium
    elif proton_number == 89:
        return 227 # Actinium
    elif proton_number == 90:
        return 232 # Thorium
    elif proton_number == 91:
        return 231 # Protactinium
    elif proton_number == 92:
        return 238 # Uranium
    else:
        sys.exit(1)        