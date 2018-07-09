#!/usr/bin/env python

"""Compute the solar magnetic field as modeled in:
Akasofu, S.-I., Gray, P., & Lee, L. 1980, Planetary Space Science, 28, 609
(1) Solar Dipole
(2) Sunspot Dipoles
(3) Solar Dynamo
(4) Ring Current
Coordinate system: (x,y,z) Sun === (0,0,0), Earth === (1,0,0)
"""

import numpy as np

import Transform

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Prototype"

###########################################################################    



###########################################################################    


def lorentzBackBoost( En_eV, mass, angle_rad ):
    """Returns the original energy of a photon measured from the reference frame of the Sun.
    That is, back-boosted from the energy encountered in the frame of a traveling nucleus En_eV 
    in [electronVolts] with mass [atomic mass units].  The relative angle between
    photon and nucleus momenta is given by angle_rad in [radians].
    Head-on, angle_rad = 0; same direction, angle_rad = pi.
    ref. 1960 Geramsimova and Zatsepin
    """
    amu2mev = 931.5 # [MeV/c**2] e.g. carbon has an atomic mass of 12, or 12 * amu2mev = 11,178 MeV/c**2
    mass_eV = mass * amu2mev * 1e6 # [eV / c**2]
    
    lorentz_gamma = En_eV / mass_eV
    
    return En_eV / ( 2. * lorentz_gamma * np.cos(angle_rad / 2.)**2 )

def solarPhotonDensity( rs_AU, Es_eV ):
    """Returns the differential solar photon number density dn/dE in 
    [number / eV * cm**3] given a radial distance from the Sun, rs_AU in
    [astronomical units] and solar photon energy Es_eV in [electronVolts] 
    as measured in the reference frame of the Sun.
    Black body spectrum with T = 5770 K.
    """
    return 7.8e7 * (1./rs_AU**2) * ( Es_eV**2 / ( np.exp(Es_eV/.5) - 1. ) )
