#!/usr/bin/env python

"""The scipy.integrate.ode module interface solves equation systems of the form, y'(t) = f(t,y).
applyForces plays the role of f(t,y), returning the computation of y'(t).
"""

import numpy as np

import Bfield

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"


def applyForces( t, Y, ratio ):
    """Computes the change in direction of relativistic-velocity (beta) for a nuclear fragment
    under the magnetic influence of the Sun as it propagates through the solar system.
    t in this context is not time.  It is the linear path displacement, s [meters], or rather l [AU]. 
    Y is the state vector (x, y, z, beta_x, beta_y, beta_z) === (position [AU], relativistic velocity [unit-less]).
    ratio === Z / E === (number protons / electronVolts) === [0..46]e-18 === [neutron..uranium]
    (i.e. intended for energy range [2..200]e18 eV, simplifying numerical assumption/requirement 
    np.dot(beta, beta) = 1. at all times).
    returns dY/dt === (beta_x, beta_y, beta_z, alpha_x, alpha_y, alpha_z) === (rel. velocity=unit-less, dv/dl=1/AU=alpha)
    """
    m_per_AU = 149597870700. # use:  position [AU] * m_per_AU = converted position [m]
    c = 299792458. # [meters / second] === speed of light    
    
    pos_x , pos_y , pos_z  = Y[0], Y[1], Y[2] # [AU]
    beta_x, beta_y, beta_z = Y[3], Y[4], Y[5] # [unit-less]
    
    pos  = np.array([ pos_x , pos_y , pos_z  ])
    beta = np.array([ beta_x, beta_y, beta_z ])
    
    ### enforce special relativity, aka numerical simplification: np.dot(beta,beta) == 1 == speed of light
    beta = beta / np.sqrt( np.dot(beta, beta) )
    beta_x, beta_y, beta_z = beta[0], beta[1], beta[2]

    def solarLorentzForce( __ratio, __pos, __beta ):
        """Computes (and returns) d/dl(beta) === (1/AU) from solar magnetic field influence.
        Derivation in the Appendix below.
        """
        B = Bfield.cartesianTesla(__pos) # [Tesla] 
        return m_per_AU * __ratio * np.cross( __beta, c*B ) # [(meters/AU) * (1/Volts) * (meters/second) * Tesla] == [1/AU]

    alpha_x, alpha_y, alpha_z = solarLorentzForce(ratio, pos, beta) # [1/AU]   
    alpha = np.array([ alpha_x, alpha_y, alpha_z ])
    
    return np.concatenate(( beta, alpha ))

"""Appendix:  Representation of the magnetic Lorentz force.

    Definitions:
       |beta| === magnitude of vector beta
        beta^ === unit vector in the direction of beta
        beta  === velocity / speed of light [unit-less] (3-vector)
        B     === applied magnetic field [Tesla] (3-vector)
        d/dt  === derivative with respect to t, [1/seconds]
        d/ds  === derivative with respect to displacement, [1 / meters]
        d/dl  === derivative with respect to displacement, [1 / AU]
        e     === charge of a proton [coulombs] (scalar constant)
        E     === energy of nuclear fragment [Joules] (scalar constant)
        EeV   === energy of nuclear fragment [electronVolts] (scalar constant)
        c     === speed of light [meters / second] (scalar constant)
        gamma === Lorentz factor === 1 / sqrt( 1 - beta**2 ) [unit-less] (scalar)
        mass  === inertial mass in rest frame [kilograms] (scalar constant)
        p     === relativistic 3-momentum of nuclear fragment [kilogram meters / second] (3-vector)
        Z     === number of protons [unit-less] (scalar constant)

    In the Sun's frame: 
        p = gamma * mass * c * beta

    Lorentz force law:
        d/dt(p) =   q   * cross-product( v, B )
        d/dt(p) = (Z*e) * cross-product( c*beta, B ) 
        
    Substitutions:
        p = (E/c) * beta
        ds = (c*|beta|) * dt

    Yields:
        (c*|beta|) d/ds( (E/c)*beta ) = (Z*e) * cross-product( c*beta, B )
        
        d/ds(beta) = (Z*e / E) * cross-product( beta^, c*B )
        
    Changing units E -> e * EeV (electronVolts carry units of Volts === Joules / coulomb)

        d/ds(beta) = (Z / EeV) * cross-product( beta^, c*B )
   
    With EeV > 2e18 electronVolts, |beta| ~= 1, therefore beta ~= beta^:
    
        d/ds(beta^) = (Z / EeV) * cross-product( beta^, c*B )
    
    Changing units ds [meters] -> dl [astronomical units] * m_per_AU:
    
        d/dl(beta^) = m_per_AU * (Z / EeV) * cross-product( beta^, c*B )
"""
