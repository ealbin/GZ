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
__status__ = "Production"


def approximateAngle(cartesian_pos):
    """Given cartesian coordinates (x, y, z) in [AU], returns the angle [radians] between
    a photon from the Sun and the Earth-bound trajectory.  This is an approximation
    as the actual angle depends on the actual trajectory which may not be directly
    Earth-bound, e.g. headed to (1, 0, 0) [AU].  But to very good approximation, this
    is true as significant deviation from Earth (aka misses the Earth) is irrelavent
    in this analysis chain. Angle is as defined in 1960 Gerasimova and Zatsepin, fig 1.   
    """
    cart_sun   = np.asarray([0, 0, 0]) # [AU]
    cart_earth = np.asarray([1, 0, 0]) # [AU]
    
    cart_ray   = cartesian_pos - cart_sun
    unit_ray   = cart_ray / np.sqrt(np.dot( cart_ray, cart_ray ))
    
    cart_traj  = cart_earth - cartesian_pos
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

def integrand(solar_e, lorentz_gamma, mass_number, sun_dist, geometry):
    """Returns the integrand of the energy integral for computing the pdf.
    See "pdf()" below.
    """    
    return ( solarPhotonDensity(sun_dist, solar_e) 
             * CrossSection.singleNucleon(mass_number, lorentz_gamma * geometry * solar_e)
             * geometry ) # [probability / centimeter * electronVolt]

def pdf(cartesian_pos, mass_number, energy_eV):
    """Returns the probility density [probability / meter] for single nucleon ejection 
    at cartesian location (x, y, z) [AU] from parent nucleide (mass_number) 
    traveling with energy (energy_eV).
    """
    amu2mev = 931.5 # [MeV/c**2] e.g. carbon has an atomic mass of 12, or 12 * amu2mev = 11,178 MeV/c**2
    mass_eV = mass_number * amu2mev * 1e6 # [eV / c**2]
    lorentz_gamma = energy_eV / mass_eV

    alpha_rad = approximateAngle(cartesian_pos)
    geometry  = 2. * np.cos(alpha_rad / 2.)**2
    
    sun_dist  = np.sqrt( np.dot( cartesian_pos, cartesian_pos ) )
    
    # analytical limits of integration are 0 to infinity [electronVolts], however the solar blackbody
    # spectrum is negligible by 10 [eV].  An upper limit of 100 [eV] is performed, as choices of a higher limit
    # will result in an undersampled integrand between 0 and 1 [eV] (the most important part). 
    pdf_cm, err = integrate.quad( integrand, 0, 100, 
                                  args=(lorentz_gamma, mass_number, sun_dist, geometry) ) 
                                  # [probability / centimeter]
        
    return pdf_cm * 100. # [probility / meter]
