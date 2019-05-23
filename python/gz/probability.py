#!/usr/bin/env python3

"""Computes the attenuation length [AU] of photodissentegration.
Coordinate system: (x,y,z) Sun === (0,0,0), Earth === (1,0,0)
"""

__project__     = 'GZ Paper'
__version__     = 'v1.0'
__objective__   = 'Phenominology'
__institution__ = 'University of California, Irvine'
__department__  = 'Physics and Astronomy'
__author__      = 'Eric Albin'
__email__       = 'Eric.K.Albin@gmail.com'
__updated__     = '13 May 2019'


import numpy as np

from scipy import integrate

from . import units
from . import cross_section
from . import photon_field


def oneOrMore(atten_length, distance):
    """Returns the probability of a process with attenuation length [AU] over
    a distance [AU].
    """
    p0 = np.exp(-distance / atten_length)
    return 1. - p0


def random(x, pdf, size, seed=None, plottables=False, CDF=False):
    x = np.asarray(x)
    pdf = np.asarray(pdf)

    cdf = np.zeros(x.size)
    for i in range(x.size):
        if (i == 0):
            continue
        cdf[i-1] = integrate.simps(pdf[:i], x=x[:i])
    cdf[-1] = cdf[-2]
    
    if (seed is not None):
        np.random.seed(seed)
    val = (cdf[-1] - cdf[0]) * np.random.random(size) + cdf[0]

    out = []
    for v in val:
        for i, c in enumerate(cdf):
            if (c < v):
                continue
            if (c == v):
                out.append(x)
                break
            c1 = cdf[i-1]
            c2 = c
            x1 = x[i-1]
            x2 = x[i]
            slope = (x2 - x1) / (c2 - c1)
            out.append(slope * (v - c1) + x1)
            break
    
    if (size == 1):
        out = out[0]
        
    if (plottables):
        if (CDF):
            out = out, x, pdf / cdf[-1], cdf
        else:
            out = out, x, pdf / cdf[-1]
    elif (CDF):
        out = out, x, cdf
        
    return out    


class Solar:    

    ENERGY_LOW_EXPONENT  = -15.
    ENERGY_HIGH_EXPONENT = 3.
    ENERGY_SAMPLES       = 2000
    
    def integrand(photon_energy, lorentz_gamma, mass_number, dist_sun, geo_factor, cross_section):
        """Returns the integrand of the energy integral for computing the attenuation.
        See "attenuation()" below.
        """    
        density = photon_field.Solar.dNdE(dist_sun, photon_energy)
        x_section = cross_section(None, lorentz_gamma * geo_factor * photon_energy, mass_number=mass_number)
        return density * x_section * geo_factor # [probability / centimeter * electronVolt]

    
    def attenuation(position, beta, proton_number, nuclide_energy, 
            cross_section=cross_section.Photodissociation.singleNucleon, 
            mass_number=None, algorithm='simps'):     
        """Returns the attenuation length [AU] for process specified by 
        cross_section parameter (single nucleon ejection by default) for
        cartesian position (x, y, z) [AU] from parent nucleide (proton_number or mass_number) 
        traveling with energy (nuclide_eV) heading in direction beta (bx, by ,bz).
        See comments below regarding 'algorithm'.  In short, 'simps' is 10x slower but accurate,
        'quad' is 10x faster but less accurate.
        """
    
        if (algorithm != 'simps' and algorithm != 'quad'):
            print('invalid algorithm: choose "simps" or "quad"')
            return
        
        if (mass_number == None):
            mass_number = units.Nuclide.mass_number(proton_number)
                
        mass_eV = mass_number * units.Change.amu_to_eV # [eV / c**2]
        lorentz_gamma = nuclide_energy / mass_eV
    
        dist_sun = np.sqrt( np.dot(position, position) )
        r_hat = position / dist_sun
        beta = beta / np.sqrt( np.dot(beta, beta) )

        dot = np.dot( -r_hat, beta)
        if (dot > 1. and np.isclose(dot, 1.)):
            dot = 1.        
        if (dot < -1. and np.isclose(dot, -1.)):
            dot = -1.            
        alpha_radians = np.arccos(dot)
        geo_factor = 2. * np.cos(alpha_radians / 2.)**2
        
        # Analytical limits of integration are 0 to infinity [electronVolts], however the solar blackbody
        # spectrum is negligible by 10 [eV].  An upper limit of 100 [eV] is performed.
        # Using an algorithm such as quad produces very similar (within ~20%) results to a sampled algorithm 
        # like simps, however I believe a well sampled simps result is closer to the true value as quad (et al)
        # tends to undersample the integrand between 0 and 1 [eV] (aka the most important part).
        # The downside is it is a few orders of magnitude slower than quad.
        if (algorithm == 'simps'):
            # good 5 digit precision at 1000 samples
            e_samples = np.logspace(Solar.ENERGY_LOW_EXPONENT, Solar.ENERGY_HIGH_EXPONENT, Solar.ENERGY_SAMPLES) 
            i_samples = np.zeros(e_samples.size)
            
            for i, e in enumerate(e_samples):
                i_samples[i] = Solar.integrand(e, lorentz_gamma, mass_number, dist_sun, geo_factor, cross_section) 
                # [probability / centimeter * electronVolt]

            atten_cm = 1. / integrate.simps( i_samples, x=e_samples ) # [centimeters]
            
        else: # algorithm == 'quad'
            # upper limit capped at 100 eV instead of infinity to avoid undersampling:
            atten_cm, err = integrate.quad( Solar.integrand, 0, 100, 
                                          args=(lorentz_gamma, mass_number, dist_sun, geo_factor, cross_section) ) 
            atten_cm = 1. / atten_cm # [centimeters]
            
        atten_m = atten_cm / 100. # [meters]
        return atten_m * units.Change.meter_to_AU # [AU]
    
    
    def get_photon(position, beta, proton_number, nuclide_energy, 
                   cross_section=cross_section.Photodissociation.singleNucleon, 
                   mass_number=None, seed=None, size=1, plottables=False, CDF=False):
        
        if (mass_number == None):
            mass_number = units.Nuclide.mass_number(proton_number)

        mass_eV = mass_number * units.Change.amu_to_eV # [eV / c**2]
        lorentz_gamma = nuclide_energy / mass_eV
    
        dist_sun = np.sqrt( np.dot(position, position) )
        r_hat = position / dist_sun
        beta = beta / np.sqrt( np.dot(beta, beta) )
    
        dot = np.dot( -r_hat, beta)
        if (dot > 1. and np.isclose(dot, 1.)):
            dot = 1.        
        if (dot < -1. and np.isclose(dot, -1.)):
            dot = -1.        
        alpha_radians = np.arccos(dot)
        geo_factor = 2. * np.cos(alpha_radians / 2.)**2
        
        e_samples = np.logspace(Solar.ENERGY_LOW_EXPONENT, Solar.ENERGY_HIGH_EXPONENT, Solar.ENERGY_SAMPLES) 
        i_samples = np.zeros(e_samples.size)
        for i, e in enumerate(e_samples):
            i_samples[i] = Solar.integrand(e, lorentz_gamma, mass_number, dist_sun, geo_factor, cross_section) 

        return random(e_samples, i_samples, size, seed=seed, plottables=plottables, CDF=CDF)