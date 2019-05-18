#!/usr/bin/env python

"""Compute the solar magnetic field as modeled in:
Akasofu, S.-I., Gray, P., & Lee, L. 1980, Planetary Space Science, 28, 609
(1) Solar Dipole
(2) Sunspot Dipoles
(3) Solar Dynamo
(4) Ring Current
Coordinate system: (x,y,z) Sun === (0,0,0), Earth === (1,0,0)
"""

import numpy as np

import Constants
import Transform

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"


###########################################################################    

### parametric constants, ref. Akasofu, Gray & Lee (1980):
Bs = 2.      # [Gauss]
Bo = -3.5e-5 # [Gauss]
Bt =  3.5e-5 # [Gauss] 
Bd = 1000.   # [Gauss]
Ro = 0.00465 # Radius of the Sun [astronomical units]
Rd = 0.1*Ro  # Sunspot dipole radius [astronomical units]
po = 1.      # [astronomical units]

###########################################################################    

def solarDipole(cartesian_pos):
    """Compute the solar dipole component of the field model given 
    cartesian position in [astronomical units].  
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    polar_pos = Transform.cartesian2polar(cartesian_pos)['rtz']
    rho       = polar_pos[0] # [astronomical units]
    theta     = polar_pos[1] # [radians]
    z         = polar_pos[2] # [astronomical units]
    
    ## B_rho [Gauss]
    B_rho = 0
    if np.abs(z) > 0:
        B_rho = -(3./2.) * (Bs * Ro**3) * rho * z * (z**2 + rho**2)**(-5./2.)    

    ## B_theta [Gauss]
    B_theta  = 0
    
    ## B_z [Gauss]
    B_z  = 0
    if np.abs(rho) > 0:
       B_z = (1./2.) * (Bs * Ro**3) * (rho**2 - 2*(z**2)) * (z**2 + rho**2)**(-5./2.)
    
    polar_B     = np.array([ B_rho, B_theta, B_z ])
    cartesian_B = Transform.polar2cartesian(polar_pos, vec=polar_B)['vec']
    return cartesian_B # [Gauss]


def solarSunspot(cartesian_pos):
    """Compute the solar sunspot component of the field model given 
    cartesian position in [astronomical units].
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    x = cartesian_pos[0] # [astronomical units]
    y = cartesian_pos[1] # [astronomical units]
    z = cartesian_pos[2] # [astronomical units]
    
    N_dipoles = 180
    dipole_thetas = np.linspace(0, 360, N_dipoles, endpoint=False) * np.pi / 180. # [radians]
    sumB_x = 0
    sumB_y = 0
    sumB_z = 0
    for dipole_theta in dipole_thetas:
        dipole_x = Rd * np.cos(dipole_theta) # [astronomical units]
        dipole_y = Rd * np.sin(dipole_theta) # [astronomical units]
        dipole_z = 0                         # [astronomical units]
        
        ## relative distance from dipole to field point
        rel_x = x - dipole_x # [astronomical units]
        rel_y = y - dipole_y # [astronomical units]
        rel_z = z - dipole_z # [astronomical units]

        rel_cartesian = np.array([ rel_x, rel_y, rel_z ])
        rel_polar     = Transform.cartesian2polar(rel_cartesian)['rtz']
        rho   = rel_polar[0] # [astronomical units]
        theta = rel_polar[1] # [radians]
        z     = rel_polar[2] # [astronomical units]
    
        ## B_rho [Gauss]
        B_rho = 0
        if np.abs(z) > 0:
            B_rho = -(3./2.) * (Bd * Rd**3) * rho * z * (z**2 + rho**2)**(-5./2.)    

        ## B_theta [Gauss]
        B_theta  = 0
    
        ## B_z [Gauss]
        B_z  = 0
        if np.abs(rho) > 0:
           B_z = (1./2.) * (Bd * Rd**3) * (rho**2 - 2*(z**2)) * (z**2 + rho**2)**(-5./2.)
            
        polar_B     = np.array([ B_rho, B_theta, B_z ])
        cartesian_B = Transform.polar2cartesian(rel_polar, vec=polar_B)['vec']

        sumB_x += cartesian_B[0]
        sumB_y += cartesian_B[1]
        sumB_z += cartesian_B[2]
        
    return np.array([ sumB_x, sumB_y, sumB_z ]) # [Gauss]    


def solarDynamo(cartesian_pos):
    """Compute the solar dynamo component of the field model given 
    cartesian position in [astronomical units].
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    polar_pos = Transform.cartesian2polar(cartesian_pos)['rtz']
    rho       = polar_pos[0] # [astronomical units]
    theta     = polar_pos[1] # [radians]
    z         = polar_pos[2] # [astronomical units]    

    ## B_rho [Gauss]
    B_rho = 0

    ## B_theta [Gauss]
    B_theta  = 0
    if np.abs(rho) > 0:
        B_theta  = (Bt * po) / float(rho)
    if z < 0:
        B_theta *= -1.
    
    ## B_z [Gauss]
    B_z = 0
    
    polar_B     = np.array([ B_rho, B_theta, B_z ])
    cartesian_B = Transform.polar2cartesian(polar_pos, vec=polar_B)['vec']
    return cartesian_B # [Gauss]    


### OPTIONAL TODO
def solarRingAGL(cartesian_pos):
    """Compute the solar ring component of the field model given 
    cartesian position in [astronomical units].
    Follows the approximation made in Akasofu, Gray & Lee (1980).
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    polar_pos = Transform.cartesian2polar(cartesian_pos)['rtz']
    rho       = polar_pos[0] # [astronomical units]
    theta     = polar_pos[1] # [radians]
    z         = polar_pos[2] # [astronomical units]    

    print("DON'T CALL ME - I'M NOT IMPLEMENTED YET")
    return np.array([ 0, 0, 0 ])
    

def solarRingEMR(cartesian_pos):
    """Compute the solar ring component of the field model given 
    cartesian position in [astronomical units].
    Follows the approximation made in Epele, Mollerach & Roulet (1999).
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    polar_pos = Transform.cartesian2polar(cartesian_pos)['rtz']
    rho       = polar_pos[0] # [astronomical units]
    theta     = polar_pos[1] # [radians]
    z         = polar_pos[2] # [astronomical units]    
    ## B_rho [Gauss]
    B_rho = 0
    if np.abs(rho) > 0:
        B_rho = (Bo * po**2) * rho * (z**2 + rho**2)**(-3./2.)
    if z < 0:
        B_rho *= -1.

    ## B_theta [Gauss]
    B_theta = 0

    ## B_z [Gauss]
    B_z = 0
    if np.abs(rho) > 0:
       B_z = (Bo * po**2) * np.abs(z) * (z**2 + rho**2)**(-3./2.)

    polar_B     = np.array([ B_rho, B_theta, B_z ])
    cartesian_B = Transform.polar2cartesian(polar_pos, vec=polar_B)['vec']
    return cartesian_B # [Gauss]    


### OPTIONAL TODO
def solarRingExact(cartesian_pos):
    """Compute the solar ring component of the field model given 
    cartesian position in [astronomical units].
    Follows the exact integral formulation in Akasofu, Gray & Lee (1980).
    returns a magnetic field density vector in cartesian coordinates in Gauss.
    """
    polar_pos = Transform.cartesian2polar(cartesian_pos)['rtz']
    rho       = polar_pos[0] # [astronomical units]
    theta     = polar_pos[1] # [radians]
    z         = polar_pos[2] # [astronomical units]    

    print("DON'T CALL ME - I'M NOT IMPLEMENTED YET")
    return np.array([ 0, 0, 0 ])    


def sumBfieldGauss(cartesian_pos):
    """Compute the total cartesian compoents of the solar magnetic field 
    given cartesian position in [astronomical units].
    Uses the EMR approximation for the solar ring field.
    returns a magnetic field density vector in Gauss.
    """
    B_dipole  = solarDipole(cartesian_pos)    # [Gauss]
    B_sunspot = solarSunspot(cartesian_pos)   # [Gauss]
    B_dynamo  = solarDynamo(cartesian_pos)    # [Gauss]
    B_ring    = solarRingEMR(cartesian_pos)   # [Gauss]
    B_total   = B_dipole + B_sunspot + B_dynamo + B_ring
    return B_total # [Gauss]
    
def sumBfieldTesla(cartesian_pos):
    """Compute the total cartesian compoents of the solar magnetic field 
    given cartesian position in [astronomical units].
    Uses the Epele approximation for the solar ring field.
    returns a magnetic field density vector in Tesla.
    """
    B_total  = sumBfieldGauss(cartesian_pos) # [Gauss]
    return Constants.Gauss2Tesla(B_total) # [Tesla]
