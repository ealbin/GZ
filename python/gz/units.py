#!/usr/bin/env python3

"""Wrapper for physical constants and conversions
"""

__project__     = 'GZ Paper'
__version__     = 'v1.0'
__objective__   = 'Phenominology'
__institution__ = 'University of California, Irvine'
__department__  = 'Physics and Astronomy'
__author__      = 'Eric Albin'
__email__       = 'Eric.K.Albin@gmail.com'
__updated__     = '13 May 2019'


class SI:
    lightspeed   = 299_792_458.  # speed of light [meters / second]
    
    radius_earth =   6_378_100.  # radius of Earth [meters]
    radius_sun   = 695_508_000.  # radius of the Sun [meters]
    

class Change:
    AU_to_meter    = 149_597_870_700.  # [meters / astronomical unit]
    meter_to_AU    = 1. / AU_to_meter
    
    amu_to_eV      = 9_314_940_954.  # [eV/c**2] mass of 1 atomic mass unit [amu] 
    eV_to_amu      = 1. / amu_to_eV
    
    tesla_to_gauss = 10_000.  # magnetic field [Gauss] in 1 [Tesla]
    gauss_to_tesla = 1. / tesla_to_gauss

    barn_to_cm2    = 1e-24  # area [barn] in cm**2
    cm2_to_barn    = 1. / barn_to_cm2


class Nuclide:

    def neutron_number(proton_number):
        """ Returns average number of neutrons for a given proton_number
        """
        A = Nuclide.mass_number(proton_number)
        if (A == None):
            return None
        return A - proton_number
    
    
    def mass_number(proton_number):
        """ Returns an integer-rounded average mass_number for a given proton_number
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
            return None        