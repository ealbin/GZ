#!/bin/env python3

import Constants
import Dynamics
import Bfield

import numpy as np
import sys


#################################################################

def GetZenith(l, b):
    """
    l = degrees longitude (0 points in +x direction)
    b = degrees latitude (0 is the x-y plane)
    returns solar system coordinates and zenith unit vector
    """
    Re = Constants.RadiusAU('earth')
    l_rad = l * np.pi / 180.
    b_rad = b * np.pi / 180.
    
    zenith_x = np.cos(b_rad) * np.cos(l_rad)
    zenith_y = np.cos(b_rad) * np.sin(l_rad)
    zenith_z = np.sin(b_rad)
    
    zenith = np.asarray([zenith_x, zenith_y, zenith_z])
    earth  = Re * zenith
    solar  = np.asarray([1,0,0]) + earth
    
    return solar, zenith

#################################################################

def GetStartPoints(l, b, n, radius=10, acceptance=90):
    """
    l = degrees longitude (0 points in the +x direction)
    b = degrees latitude (0 is the x-y plane)
    n = number of points to get
    radius = distance away from sun in AU
    acceptance = degrees off of zenith, should be 90 (earth surface plane) or less
    returns start location and initial heading vectors
    """
    def UnitVectors(n):
        """
        lr = radians
        br = radians
        return unit vector
        """
        lr = 2.*np.pi * np.random.random(size=n)
        br = np.pi * np.random.random(size=n) - np.pi/2.
        x = np.cos(br) * np.cos(lr)
        y = np.cos(br) * np.sin(lr)
        z = np.sin(br)
        return np.column_stack((x,y,z))
    
    earth, zenith = GetZenith(l, b)
    
    points = np.array([]).reshape(0,3)
    trajectory = np.array([]).reshape(0,3)
    
    while (len(points) < n):
        p = radius * UnitVectors(n - len(points))
        t = earth - p
            
        # make unit vector
        t = t / np.sqrt( (t*t).sum(1) ).reshape(1,-1).T
            
        angles = 180. / np.pi * np.arccos( np.dot(t, -zenith) )
    
        # keep only points forward of the earth surface plane
        condition = (angles <= acceptance)
    
        points = np.concatenate((points, p[condition]))
        trajectory = np.concatenate((trajectory, t[condition]))
    
    return points, trajectory

#################################################################

def HaversineNormalized(lat1, lon1, lat2, lon2):
    """ Normalized the earth radius aka 1 = Re
    """
    lat1 *= np.pi / 180.
    lon1 *= np.pi / 180.
    lat2 *= np.pi / 180.
    lon2 *= np.pi / 180.
    
    part1 = np.sin( (lat2 - lat1) / 2. )**2.
    part2 = np.cos(lat1) * np.cos(lat2)
    part3 = np.sin( (lon2 - lon1) / 2. )**2.

    return 2. * np.arcsin( np.sqrt(part1 + part2 * part3) )

#################################################################

def Trajectories(l, b, ratio, n=1, radius=10, acceptance=90):
    
    m_per_AU = Constants.AU2meters(1.) # meters in 1 AU
    c = Constants.lightspeed_m_per_s() # speed of light in [meters / second]
    Re = Constants.RadiusAU('earth')
    E = np.asarray([1,0,0])
     
    out = []
    
    points, heading = GetStartPoints(l, b, n, radius=radius, acceptance=acceptance)
    print('Starting points (x,y,z) [AU]:')
    print(points)
    for point, heading in zip(points, heading):
        
        position = np.asarray(point)
        beta = np.asarray(heading)
        
        r_earth = position - E
        dist_earth = np.sqrt( np.dot(r_earth, r_earth) )
        limit = dist_earth + 2.*Re
        
        r_hat = r_earth / dist_earth
        cosE = np.dot(beta, -r_hat)
        x = (dist_earth * cosE) - np.sqrt(Re**2 - dist_earth**2 * (1. - cosE*cosE))
        R = r_earth + (x * beta)
        f = np.sqrt(np.dot(R,R)) / Re
        u = R / np.sqrt(np.dot(R,R))
        ll = (180./np.pi) * np.arctan2(u[1], u[0])
        bb = (180./np.pi) * np.arctan2(u[2], np.sqrt(u[0]*u[0] + u[1]*u[1]))
        
        step_max = limit / 1e4 #np.log10(limit / 1e4)
        step_mid = 5. * Re 
        
        near_earth = 5. * Re
        step_slope1 = (step_max - step_mid) / limit        
        step = step_slope1 * dist_earth + step_mid #np.power(10., step_slope * dist_earth + step_min)
        
        dist = 0.
        while (dist < limit):
    
            beta = beta / np.sqrt( np.dot(beta, beta) )
                        
            B = Bfield.cartesianTesla(position) # [Tesla] 
            alpha = m_per_AU * ratio * np.cross( beta, c*B )
    
            beta = beta + (alpha * step)
            position = position + (beta * step)
            
            r_earth = position - E
            dist_earth = np.sqrt( np.dot(r_earth, r_earth) )
            if (dist_earth > near_earth):
                step = step_slope1 * dist_earth + step_mid
            else:
                r_hat = r_earth / dist_earth
                cosE = np.dot(beta, -r_hat)
                ratical = Re**2 - dist_earth**2 * (1. - cosE*cosE)
                if (ratical < 0.):
                    #print('misses the earth, ending this simulation')
                    break
                x = (dist_earth * cosE) - np.sqrt(ratical)
                R = r_earth + (x * beta)
                f = np.sqrt(np.dot(R,R)) / Re
                if (not np.isclose(f, 1.)):
                    #print('extrapolation of trajectory failed, skipping')
                    break
                u = R / np.sqrt(np.dot(R,R))
                ll = (180./np.pi) * np.arctan2(u[1], u[0])
                bb = (180./np.pi) * np.arctan2(u[2], np.sqrt(u[0]*u[0] + u[1]*u[1]))
                out.append(HaversineNormalized(b, l, bb, ll) * Constants.RadiusMeters('earth'))
                #print(' finished', end='', flush=True)
                break
            
            dist += step
            
        #print()
    
    print('Distances in meters:')
    return out

if (len(sys.argv) < 5 or len(sys.argv) > 8):
	print('usage: longitude(deg), latitude(deg), Z, E, optional (n, radius, acceptance)')
else:
	seed = np.random.randint(0, 999999999)
	print('Using random seed: ' + str(seed))
	np.random.seed(seed)

	l = float(sys.argv[1])
	b = float(sys.argv[2])
	Z = float(sys.argv[3])
	E = float(sys.argv[4])

	print('Target longitude(l): ' + str(l) + ' [deg]')
	print('Target latitude(b): ' + str(b) + ' [deg]')
	print('Atomic number(Z): ' + str(Z))
	print('Energy(E): ' + str(E) + ' [eV]')

	ratio = Z / E
	n = 1 
	radius = 10 # AU
	acceptance = 90 # deg
	
	if (len(sys.argv) > 5):
		n = int(sys.argv[5])
	if (len(sys.argv) > 6):
		radius = float(sys.argv[6])
	if (len(sys.argv) > 7):
		acceptance = float(sys.argv[7])

	print('Radius(R): ' + str(radius) + ' [AU]')
	print('Acceptance(A): ' + str(acceptance) + ' [deg]')

	print( Trajectories(l, b, ratio, n=n, radius=radius, acceptance=acceptance) )