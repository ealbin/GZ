#!/usr/bin/env python3

"""Special relativity
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

from . import units

def gamma(energy, mass):
    """ Energy in eV, mass in eV/c*2
    """
    return energy / mass

def beta(gamma):
    return np.sqrt(1. - 1. / gamma)

def momentum_mag(energy, mass):
    return np.sqrt(energy*energy - mass*mass)

def momentum(energy, mass, direction):
    """ direction is a unit vector
    """
    return momentum_mag(energy, mass) * np.asarray(direction)

def theta(e1, p1, e2, p2, e3, m3, e4, m4):
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    
    e12 = (e1 + e2)**2
    p12 = np.dot(p1 + p2, p1 + p2)
    
    e34 = (e3 + e4)**2
    p3_mag = momentum_mag(e3, m3)
    p4_mag = momentum_mag(e4, m4)
    
    cosTheta = ( e34 - (e12 - p12) - p3_mag**2 - p4_mag**2 ) / (2. * p3_mag * p4_mag)
    if (cosTheta > 1. and np.isclose(cosTheta, 1.)):
        cosTheta = 1.
    if (cosTheta < -1. and np.isclose(cosTheta, -1.)):
        cosTheta = -1.
    return np.arccos(cosTheta)