#!/usr/bin/env python3

"""Computes the photon field density in [number / (eV * cm**3)]
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

from . import coordinates
from . import units

class Solar:
    
    def earthShadow(position):
        Re = units.SI.radius_earth * units.Change.meter_to_AU
        Rs = units.SI.radius_sun * units.Change.meter_to_AU
        
        earth = coordinates.Cartesian.earth
        sun = coordinates.Cartesian.sun
        
        p2earth = earth - position
        p2sun = sun - position
        
        p2earth_dist = np.sqrt(np.dot(p2earth, p2earth))
        p2sun_dist = np.sqrt(np.dot(p2sun, p2sun))
        
        # Inside Earth
        if (p2earth_dist < Re):
            return 0.
        
        # On the darkside of the Earth
        if (p2earth_dist <= Re and position[0] > earth[0]):
            return 0.
        
        # Earth is behind the Sun
        if (p2earth_dist > p2sun_dist):
            return 1.
        
        earth_sun_angle = np.arccos( np.dot(p2earth, p2sun) / ( p2earth_dist * p2sun_dist ) )
        earth_angle = np.arcsin(Re / p2earth_dist)
        sun_angle = np.arcsin(Rs / p2sun_dist)
        
        # Apparent Earth radius
        re = Re
        
        # Apparent Sun radius
        rs = p2earth_dist * np.sin(sun_angle)
        
        # Apparent distance between Earth and Sun objects
        d = p2earth_dist * np.sqrt(2. * (1. - np.cos(earth_sun_angle)) )
        
        # Earth is not obscuring the Sun
        if (re <= d - rs):
            return 1.
        
        # Earth and Sun perfectly aligned
        if (earth_sun_angle == 0.):
            if (re > rs):
                return 0.
            else:
                return 1. - (re*re)/(rs*rs)
        
        # Earth fully inside Sun, or vise-versa
        if (rs > re):
            rbig = rs
            rsmall = re
        else:
            rbig = re
            rsmall = rs        
        if (d + rsmall <= rbig):
            if (re > rs):
                return 0.
            else:
                return 1. - (re*re)/(rs*rs)
        
        # Area overlapping
        def arg(r1, r2):
            out = (d*d + r1*r1 - r2*r2) / (2. * d * r1)
            return out
        
        a1 = re*re * np.arccos(arg(re, rs))
        a2 = rs*rs * np.arccos(arg(rs, re))
        a3 = (-d + re + rs) * (d + re - rs) * (d - re + rs) * (d + re + rs)
        a3 = .5 * np.sqrt(a3)
        A = a1 + a2 - a3
                      
        # Fraction of Sun showing
        Asun = np.pi * rs*rs
        return 1. - (A / Asun)
    
    
    def dNdE(distance_AU, energy_eV, position=None):
        """Returns the differential solar photon number density dn/dE in 
        [number / eV * cm**3] given a radial distance from the Sun, distance_AU in
        [astronomical units] and solar photon energy energy_eV in [electronVolts] 
        as measured in the reference frame of the Sun.
        Black body spectrum with T = 5770 K.
        """    
        if (energy_eV == 0.):
            return 0.
        scale = 7.8e7  
        r_dependence = 1. / distance_AU**2
        exponent = energy_eV / .5
        if (np.abs(exponent) > 100.):
            return 0.
        e_dependence = energy_eV**2 / ( np.exp(exponent) - 1. )
        
        shadow = 1.
        if (position is not None):
            shadow = Solar.earthShadow(position)
        return shadow * scale * r_dependence * e_dependence


class CMB:
    
    # TODO:
    def dNdE(energy_eV):
        if (energy_eV == 0.):
            return 0.
        scale = 1. / np.pi**2 # / (hbar c)**3
        kT = 1. # = kB * 2.725 K
        exponent = energy_eV / kT
        if (np.abs(exponent) > 100.):
            return 0.
        e_dependence = ( energy_eV**2 / np.exp(exponent) - 1. )
        return scale * e_dependence