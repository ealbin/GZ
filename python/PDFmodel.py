#!/usr/bin/env python

"""Computes the probability density [per meter] of photodissentegration.
Coordinate system: (x,y,z) Sun === (0,0,0), Earth === (1,0,0)
"""

import numpy as np
import sys

from scipy import integrate

import Constants
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

def pdf(cartesian_pos, mass_number, energy_eV, algorithm='simps'):
    """Returns the probility density [probability / meter] for single nucleon ejection 
    at cartesian location (x, y, z) [AU] from parent nucleide (mass_number) 
    traveling with energy (energy_eV).
    See comments below regarding 'algorithm'.  In short, 'simps' is very slow but accurate,
    'quad' is very fast but less accurate.
    """
    mass_eV = Constants.amu2eV(mass_number) # [eV / c**2]
    lorentz_gamma = energy_eV / mass_eV

    alpha_rad = approximateAngle(cartesian_pos)
    geometry  = 2. * np.cos(alpha_rad / 2.)**2
    
    sun_dist  = np.sqrt( np.dot( cartesian_pos, cartesian_pos ) )
    
    # Analytical limits of integration are 0 to infinity [electronVolts], however the solar blackbody
    # spectrum is negligible by 10 [eV].  An upper limit of 100 [eV] is performed.
    # Using an algorithm such as quad produces very similar (within ~20%) results to a sampled algorithm 
    # like simps, however I believe a well sampled simps result is closer to the true value as quad (et al)
    # tends to undersample the integrand between 0 and 1 [eV] (aka the most important part).
    # The downside is it is a few orders of magnitude slower than quad.
    if algorithm == 'simps':
        e_samples = np.logspace(-4, 2, 1000) # good 5 digit precision at 1000 samples
        p_samples = np.zeros(e_samples.size)
        for i, e in enumerate(e_samples):
            p_samples[i] = integrand(e, lorentz_gamma, mass_number, sun_dist, geometry) 
            # [probability / centimeter * electronVolt]
        pdf_cm = integrate.simps( p_samples, x=e_samples ) # [proability / centimeter]
        
    elif algorithm == 'quad':
        # upper limit capped at 100 eV instead of infinity to avoid undersampling:
        pdf_cm, err = integrate.quad( integrand, 0, 100, 
                                      args=(lorentz_gamma, mass_number, sun_dist, geometry) ) 
                                      # [probability / centimeter]
    else:
        print 'wrong algorithm: choose "simps" or "quad"'
        sys.exit(1)
        
    return pdf_cm * 100. # [probility / meter]
