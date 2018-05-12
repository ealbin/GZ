import SolarMagneticModel
import numpy as np
from scipy.integrate import ode
"""scipy.integrate.ode module interface solves equation systems of the form, y'(t) = f(t,y)

applyForces plays the role of f(t,y), returning the computation of y'(t)
"""

def applyForces( t, Y, ratio ):
    """Computes the change in direction of relativistic-velocity (beta) for a nuclear fragment
    under the magnetic influence of the Sun as it propagates through the solar system.
    t in this context is not time.  It is the linear path displacement, s. 
    Y is the state vector (x, y, z, beta_x, beta_y, beta_z) === (position, velocity).
    ratio === Z / E === (number protons / electronVolts) === [0..46]e-18 === [neutron..uranium]
    (i.e. intended for energy range [2..200]e18 eV, simplifying numerical assumption/requirement 
    np.dot(beta, beta) = 1. at all times).
    returns dY/dt === (beta_x, beta_y, beta_z, freq_x, freq_y, freq_z) === (velocity, dv/ds=1/time=freq)
    """
    pos_x , pos_y , pos_z  = Y[0], Y[1], Y[2]
    beta_x, beta_y, beta_z = Y[3], Y[4], Y[5]
    
    pos  = np.array([ pos_x , pos_y , pos_z  ])
    beta = np.array([ beta_x, beta_y, beta_z ])

    ### enforce special relativity, aka numerical simplification: np.dot(beta,beta) == 1 == speed of light
    beta /= np.sqrt( np.dot(beta, beta) )
    beta_x, beta_y, beta_z = beta[0], beta[1], beta[2]

    """Forces
    """
    def testLorentzForce( ignored, __pos, __beta ):
        """Computes (and returns) d/ds(beta) for the uniform test field.
        rule of thumb:  r = 2.23e-20 * E / ( Z * B )
        r === radius [astronomical units]
        E === energy [electronVolts]
        Z === number of protons [unit-less]
        B === magnetic field [Tesla]
        If everything works correctly, result will be a circular orbit with
        radius ~0.0005 [AU]
        """
        c = 2.99792458e8 # [meters / second] === speed of light
        B = np.array([ 0, 0, 1 ]) # test field [Tesla]
        __ratio = 46e-18 # [1/Volts] (uranium, Z=92, E=2e18 eV)
        return __ratio * np.cross( __beta, c*B ) # [(1/Volts) * (meters/second) * Tesla] == [1/meters]

    def solarLorentzForce( __ratio, __pos, __beta ):
        """Computes (and returns) d/ds(beta) from solar magnetic field influence.
        Derivation in Appendix A below.
        """
        c = 2.99792458e8 # [meters / second] === speed of light
        B = SolarMagneticModel.sumBfieldTesla(__pos) # [Tesla] 
        return __ratio * np.cross( __beta, c*B ) # [(1/Volts) * (meters/second) * Tesla] == [1/meters]

    freq_x, freq_y, freq_z = solarLorentzForce(ratio, pos, beta)
    freq = np.array([ freq_x, freq_y, freq_z ])
    
    return np.concatenate(( beta, freq ))

"""Appendix A.
Representation of the magnetic Lorentz force.

    Definitions:
       |beta| === magnitude of vector beta
        beta^ === unit vector in the direction of beta
        beta  === velocity / speed of light [unit-less] (3-vector)
        B     === applied magnetic field [Tesla] (3-vector)
        d/dt  === derivative with respect to t
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
"""

