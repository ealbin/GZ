#!/usr/bin/env python

"""Computes the probability density [per meter] of photodissentegration.
Coordinate system: (x,y,z) Sun === (0,0,0), Earth === (1,0,0)
"""

import numpy as np

from scipy import integrate

import CrossSection

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Development"


def approximateAngle(x, y, z):
    """Given cartesian coordinates (x, y, z) in [AU], returns the angle [radians] between
    a photon from the Sun and the Earth-bound trajectory.  This is an approximation
    as the actual angle depends on the actual trajectory which may not be directly
    Earth-bound, e.g. headed to (1, 0, 0) [AU].  But to very good approximation, this
    is true as significant deviation from Earth (aka misses the Earth) is irrelavent
    in this analysis chain. Angle is as defined in 1960 Gerasimova and Zatsepin, fig 1.   
    """
    cart_sun   = np.asarray([0, 0, 0]) # [AU]
    cart_earth = np.asarray([1, 0, 0]) # [AU]
    
    cart_ray   = np.asarray([x, y, z]) - cart_sun
    unit_ray   = cart_ray / np.sqrt(np.dot( cart_ray, cart_ray ))
    
    cart_traj  = cart_earth - np.asarray([x, y, z])
    unit_traj  = cart_traj / np.sqrt(np.dot( cart_traj, cart_traj ))

    return np.arccos( np.dot( -unit_ray, unit_traj ) )

def solarPhotonDensity( rs_AU, Es_eV ):
    """Returns the differential solar photon number density dn/dE in 
    [number / eV * cm**3] given a radial distance from the Sun, rs_AU in
    [astronomical units] and solar photon energy Es_eV in [electronVolts] 
    as measured in the reference frame of the Sun.
    Black body spectrum with T = 5770 K.
    """
    return 7.8e7 * (1./rs_AU**2) * ( Es_eV**2 / ( np.exp(Es_eV/.5) - 1. ) )

def pdf(x, y, z, mass_number, energy_eV):
    """Returns the probility density [probability / meter] for single nucleon ejection 
    at cartesian location (x, y, z) [AU] from parent nucleide (mass_number) 
    traveling with energy (energy_eV).
    """
    amu2mev = 931.5 # [MeV/c**2] e.g. carbon has an atomic mass of 12, or 12 * amu2mev = 11,178 MeV/c**2
    mass_eV = mass_number * amu2mev * 1e6 # [eV / c**2]
    lorentz_gamma = energy_eV / mass_eV

    alpha_rad = approximateAngle(x, y, z)
    geometry  = 2. * np.cos(alpha_rad / 2.)**2
    
    cart_pos  = np.asarray([x, y, z])
    sun_dist  = np.sqrt( np.dot( cart_pos, cart_pos ) )
    
    pdf_cm, err = integrate.quad(lambda solar_e: solarPhotonDensity(sun_dist, solar_e) 
                                 * CrossSection.singleNucleon(mass_number, lorentz_gamma * geometry * solar_e)
                                 * geometry, 0, np.inf ) # [probability / centimeter]
    
    return pdf_cm * 100. # [probility / meter]
