#!/usr/bin/env python3

"""Transformations between coordinate systems.
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


def cartesian2polar( xyz, vec=np.array([0,0,0]) ):
    """Transform from cartesian x-y-z coordinates to polar rho-theta-z.
    Optionally also transform a cartesian vector into a polar one.
    returns dictionary 'rtz':(position) and 'vec':(transformed vector)
    """
    x = xyz[0] # [distance units]
    y = xyz[1] # [distance units]
    z = xyz[2] # [distance units]

    vec_x = vec[0] # [any unit]
    vec_y = vec[1] # [any unit]
    vec_z = vec[2] # [any unit]
    
    ### convert to polar
    rho   = np.sqrt( x**2 + y**2 ) # [distance units]
    theta = np.arctan2(y, x)       # [radians]
    z     = z                      # [distance units]

    vec_rho   =  vec_x * np.cos(theta) + vec_y * np.sin(theta)
    vec_theta = -vec_x * np.sin(theta) + vec_y * np.cos(theta)
    vec_z     =  vec_z

    return { 'rtz':np.array([ rho, theta, z ]), 
             'vec':np.array([ vec_rho, vec_theta, vec_z ]) }

def polar2cartesian( rtz, vec=np.array([0,0,0]) ):
    """Transform from polar rho-theta-z coordinates to cartesian x-y-z.
    Optionally also transform a polar vector into a cartesian one.
    returns dictionary 'xyz':(position) and 'vec':(transformed vector)
    """
    rho   = rtz[0] # [distance units]
    theta = rtz[1] # [radians]
    z     = rtz[2] # [distance units]
    
    vec_rho   = vec[0] # [any unit]
    vec_theta = vec[1] # [any unit]
    vec_z     = vec[2] # [any unit]
    
    ### convert to cartesian
    x = rho * np.cos(theta) # [distance units]
    y = rho * np.sin(theta) # [distance units]
    z = z                   # [distance units]

    vec_x = vec_rho * np.cos(theta) - vec_theta * np.sin(theta)
    vec_y = vec_rho * np.sin(theta) + vec_theta * np.cos(theta)
    vec_z = vec_z

    return { 'xyz':np.array([ x, y, z ]), 
             'vec':np.array([ vec_x, vec_y, vec_z ]) }
        
class Cartesian:
    sun   = np.asarray([0,0,0]) # [AU, AU, AU]
    earth = np.asarray([1,0,0]) # [AU, AU, AU]

class Polar:
    sun   = cartesian2polar(Cartesian.sun)['rtz']   # [AU, radian, AU]
    earth = cartesian2polar(Cartesian.earth)['rtz'] # [AU, radian, AU]

class Spherical:
    def rotate(vector, theta, phi):
        vector = np.asarray(vector, dtype=np.float64)
        v = np.sqrt(np.dot(vector, vector))
       
        u = vector / v
        u_theta = np.arccos(u[2])
        u_phi = np.arctan2(u[1], u[0])
                
        x = v * np.sin(u_theta + theta) * np.cos(u_phi + phi)
        y = v * np.sin(u_theta + theta) * np.sin(u_phi + phi)
        z = v * np.cos(u_theta + theta)
        return np.asarray([x, y, z])