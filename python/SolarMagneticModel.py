### TODO:  Bz is wrong(?) but the equations match the papers...

# intended use:  Bfield(position) (positions in AU)
# returns field strength [ Brho, Btheta, Bz ] in Tesla

import numpy as np

# magnetic units:  1 T = 1e4 G
Gauss2Tesla = 1e-4 # B[G] * Gauss2Tesla = B[T]

Bs = 2.      # G
Bo = -3.5e-5 # G
Bt =  3.5e-5 # G 

# distance units:  1 AU = 149597870700 m
Ro = 0.00465 # AU
po = 1.      # AU

def B_rho(position):
    rho   = position[0]
    theta = position[1]
    z     = position[2]
    
    dipole = 0
    if np.abs(rho) > 0 and np.abs(z) > 0:
        dipole = -(3./2.) * (Bs * Ro**3) * rho * z * (z**2 + rho**2)**(-5./2.)

    dynamo = 0

    ring   = 0
    if np.abs(rho) > 0:
        ring = (Bo * po**2) * rho * (z**2 + rho**2)**(-3./2.)
    if z < 0:
        ring *= -1.

    sunspot = 0
    return ( dipole + dynamo + ring + sunspot ) * Gauss2Tesla
        
def B_theta(position):
    rho   = position[0]
    theta = position[1]
    z     = position[2]
    
    dipole  = 0
    
    dynamo  = 0
    if np.abs(rho) > 0:
        dynamo  = (Bt * po) / rho
    if z < 0:
        dynamo *= -1.
        
    ring    = 0
    
    sunspot = 0
    return ( dipole + dynamo + ring + sunspot ) * Gauss2Tesla
    
def B_z(position):
    rho   = position[0]
    theta = position[1]
    z     = position[2]
    
    dipole  = 0
    if np.abs(rho) > 0:
       dipole  = (1./2.) * (Bs * Ro**3) * (rho**2 - 2*(z**2)) * (z**2 + rho**2)**(-5./2.)

    dynamo  = 0

    ring    = 0
    if np.abs(rho) > 0:
       ring = (Bo * po**2) * np.abs(z) * (z**2 + rho**2)**(-3./2.)

    sunspot = 0
    return ( dipole + dynamo + ring + sunspot ) * Gauss2Tesla

def Bfield(position):
    return np.array( [ B_rho(position), B_theta(position), B_z(position) ] )
