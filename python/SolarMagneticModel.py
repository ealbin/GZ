### TODO:URGENT:  solarB_z() is wrong(??) but the equations match the papers...?

"""Compute the solar magnetic field as modeled in:
Akasofu, S.-I., Gray, P., & Lee, L. 1980, Planetary Space Science, 28, 609
"""
import numpy as np

### magnetic units:  1 Telsa = 1e4 Gauss
Gauss2Tesla = 1e-4 # e.g. B-field [Gauss] * Gauss2Tesla = B-field [Tesla]

### parametric constants, ref. Akasofu, Gray & Lee (1980):
Bs = 2.      # [Gauss]
Bo = -3.5e-5 # [Gauss]
Bt =  3.5e-5 # [Gauss] 
Ro = 0.00465 # Radius of the Sun [astronomical units]
po = 1.      # [astronomical units]

###########################################################################    

def solarB_rho(rho, theta, z):
    """Compute the polar rho-component of the solar magnetic field given 
    polar position in [astronomical units].
    returns a scalar magnetic field value in Gauss.
    """
    
    ## Solar Dipole Component [Gauss]
    dipole = 0
    if np.abs(z) > 0:
        dipole = -(3./2.) * (Bs * Ro**3) * rho * z * (z**2 + rho**2)**(-5./2.)

    ## Solar Dynamo Component [Gauss]
    dynamo = 0

    ## Solar Ring Component [Gauss]
    ring   = 0
    if np.abs(rho) > 0:
        ring = (Bo * po**2) * rho * (z**2 + rho**2)**(-3./2.)
    if z < 0:
        ring *= -1.

    ## Solar Sunspot Component [Gauss]
    sunspot = 0 
    
    return ( dipole + dynamo + ring + sunspot ) # [Gauss]
        
def solarB_theta(rho, theta, z):
    """Compute the polar theta-component of the solar magnetic field given polar
    position in [astronomical units].
    returns a scalar magnetic field value in Gauss.
    """
    ## Solar Dipole Component [Gauss]
    dipole  = 0
    
    ## Solar Dynamo Component [Gauss]
    dynamo  = 0
    if np.abs(rho) > 0:
        dynamo  = (Bt * po) / rho
    if z < 0:
        dynamo *= -1.

    ## Solar Ring Component [Gauss]        
    ring    = 0
    
    ## Solar Sunspot Component [Gauss]
    sunspot = 0
    
    return ( dipole + dynamo + ring + sunspot ) # [Gauss]
    
def solarB_z(rho, theta, z):
    """Compute the polar z-component of the solar magnetic field given polar
    position in [astronomical units].
    returns a scalar magnetic field value in Gauss.
    """
    ## Solar Dipole Component [Gauss]
    dipole  = 0
    if np.abs(rho) > 0:
       dipole  = (1./2.) * (Bs * Ro**3) * (rho**2 - 2*(z**2)) * (z**2 + rho**2)**(-5./2.)

    ## Solar Dynamo Component [Gauss]
    dynamo  = 0

    ## Solar Ring Component [Gauss]
    ring    = 0
    if np.abs(rho) > 0:
       ring = (Bo * po**2) * np.abs(z) * (z**2 + rho**2)**(-3./2.)

    ## Solar Sunspot Component [Gauss]
    sunspot = 0

    return ( dipole + dynamo + ring + sunspot ) # [Gauss]

def solarBfield(position):
    """Compute the **cartesian compoents** of the solar magnetic field given **cartesian**
    position in [astronomical units].
    returns a scalar magnetic field value in **Tesla**.
    """
    x = position[0] # [astronomical units]
    y = position[1] # [astronomical units]
    z = position[2] # [astronomical units]
    
    ### convert to polar
    rho   = np.sqrt( x**2 + y**2 ) # [astronomical units]
    theta = np.arctan2(y, x)       # [radians]
    z     = z                      # [astronomical units]

    B_rho   = solarB_rho(rho, theta, z)   # [Gauss]
    B_theta = solarB_theta(rho, theta, z) # [Gauss]
    B_z     = solarB_z(rho, theta, z)     # [Gauss]

    ### convert back to cartesian
    B_x = B_rho * np.cos(theta) - B_theta * np.sin(theta)
    B_y = B_rho * np.sin(theta) + B_theta * np.cos(theta)
    B_z = B_z

    B = np.array([ B_rho, B_theta, B_z ]) # [Gauss]
    B *= Gauss2Tesla # [Tesla]

    return B # [Tesla]
