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

class Solar:
    
    def dNdE(distance_AU, energy_eV):
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
        if (exponent > 100.):
            return 0.
        e_dependence = energy_eV**2 / np.exp(exponent) - 1.
        return scale * r_dependence * e_dependence

class CMB:
    
    # TODO:
    def dNdE(energy_eV):
        if (energy_eV == 0.):
            return 0.
        scale = 1. / np.pi**2 # / (hbar c)**3
        kT = 1. # = kB * 2.725 K
        exponent = energy_eV / kT
        if (exponent > 100.):
            return 0.
        e_dependence = energy_eV**2 / np.exp(exponent) - 1.
        return scale * e_dependence