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
        scale = 7.8e7  
        r_dependence = 1. / distance_AU**2
        e_dependence = energy_eV**2 / ( np.exp(energy_eV / .5) - 1. ) 
        return scale * r_dependence * e_dependence
