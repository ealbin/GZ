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
    Y is the state vector (rho, theta, z, beta_rho, beta_theta, beta_z) === (position, velocity).
    ratio === Z / E === (number protons / electronVolts) === [0..46]e-18 === [neutron..uranium]
    (i.e. intended for energy range [2..200]e18 eV, simplifying numerical assumption/requirement 
    np.dot(beta, beta) = 1. at all times).
    returns dY/dt === (beta_rho, beta_theta, beta_z, freq_rho, freq_theta, freq_z) === (velocity, dv/ds=1/time=freq)
    """
    pos_rho , pos_theta , pos_z  = Y[0], Y[1], Y[2]
    beta_rho, beta_theta, beta_z = Y[3], Y[4], Y[5]
    
    pos  = np.array([ pos_rho , pos_theta , pos_z  ])
    beta = np.array([ beta_rho, beta_theta, beta_z ])

    ### enforce special relativity, aka numerical simplification: np.dot(beta,beta) == 1 == speed of light
    beta /= np.sqrt( np.dot(beta, beta) )
    beta_rho, beta_theta, beta_z = beta[0], beta[1], beta[2]

    """Forces
    """
    def solarLorentzForce( __ratio, __pos, __beta ):
        """Computes (and returns) dBeta_dS from solar magnetic field influence.
        Derivation in Appendix A below.
        """
        c = 3.e8 # meters / second === speed of light
        B = SolarMagneticModel.Bfield(__pos) # Tesla 
        return __ratio * np.cross( __beta, c*B ) # (1/electronVolts) * (meters/second) * Tesla == 

    freq_rho, freq_theta, freq_z = lorentzForce(ratio, pos, beta)
    freq = np.array([ freq_rho, freq_theta, freq_z ])
    
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



initial_s    = 0
initial_pos  = np.array([ 1, 0, 0 ])
initial_beta = np.array([ -.707, -.707, 0 ])
initial_conditions = np.concatenate(( initial_pos, initial_beta ))
ratio = 46 * 10**(-18)

r = ode(applyForces).set_integrator('dori5')
r.set_initial_value(initial_conditions, initial_s).set_f_params(ratio)

AU2m = 149597870700. # use:  number [AU] * AU2m = converted number [m]
final_s = 5 * AU2m   # track particle for 5 AU (in meters)
ds = 1e-3 * AU2m     # numerical stepsize, ds distance (in meters)
positions = []
while r.successful() and r.t < final_s:
    r.integrate(r.t + ds)
    positions.append(r.y[:3] / AU2m)

positions = np.array(positions)

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot3D(positions[:,0], positions[:,1]*AU2m, positions[:,2]*AU2m)

plt.show()
#plt.plot( positions[:,1], positions[:,2] )
#plt.show()
