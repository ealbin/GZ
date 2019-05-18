#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description
"""

import time

import numpy as np

import Constants
import Bfield
import PDFmodel

m_per_AU = Constants.AU2meters(1.) # meters in 1 AU
c = Constants.lightspeed_m_per_s() # speed of light in [meters / second]
Re = Constants.RadiusAU('earth') # earth radius in AU
Earth = np.asarray([1,0,0]) # earth location in AU
     


class Path:
    
    def __init__(self, position, beta, Z, E):
        """
        position: np.array(x,y,z) start position
        beta: np.array(bx,by,bz) start beta (direction of propagation)
        Z: atomic number
        E: energy in eV
        """
        self.position = position
        self.beta = beta
        self.Z = Z
        self.E = E
        self.ratio = Z / E
        
        self.start_position = position
        self.start_beta = beta
        self.distance = 0.
        
        self.__set_dist_earth()
        self.__set_dist_sun()
    
    def __set_stepsize(self):
        # TODO: default stepsize
        return 1.
    
    def __set_dist_earth(self):
        global Earth
        r = self.position - Earth
        self.dist_earth = np.sqrt(np.dot(r, r))
        
    def __set_dist_sun(self):
        r = self.position
        self.dist_sun = np.sqrt(np.dot(r, r))
    
    def propogate(self, B_override=None, step_override=None):
        """
        Propagates one step
        B_override: use this B instead of Bfield
        step_override: use this step instead of step()
        """
        global m_per_AU, c
        
        if (B_override != None):
            B = B_override
        else:
            B = Bfield.cartesianTesla(self.position) 
        
        if (step_override != None):
            step = step_override
        else:
            step = self.stepsize()
        
        alpha = m_per_AU * self.ratio * np.cross( self.beta, c*B )

        self.beta = self.beta + (alpha * step)
        self.position = self.position + (self.beta * step)

        self.distance += step
        self.__set_dist_earth()
        self.__set_dist_sun()
        
    
class Outgoing(Path):
    
    save_step = 100
    
    def __init__(self, position, beta, Z, A, E, R_limit=10.):
        """
        position: np.array(x,y,z) AU start position on earth
        beta: np.array(bx,by,bz) start beta (away from earth)
        Z: atomic number
        A: atomic mass
        E: energy in eV
        R_limit: radius [AU] of maximum propagation
        """
        global Re
        
        super().__init__(position, beta, Z, E)
        self.A = A
        self.R_limit = R_limit
        
        step_min = Re
        step_max = R_limit / 1e3
        slope = (step_max - step_min) / R_limit
        self.step_func = lambda d: slope * d + step_min
        
        self.telementry = []
        self.save_count = 0
        
    def __set_stepsize(self):
        return self.step_func(self.dist_sun)
    
    def __add_telemetry(self):
        if (self.dist_sun < 1.5 or self.save_count >= Outgoing.save_step):
            self.telementry.append([self.distance, self.pdf])     
            self.save_count = 0
        else:
            self.save_count += 1
    
    def propogate(self, B_override=None, step_override=None):
        """
        Propagates one step
        B_override: use this B instead of Bfield
        step_override: use this step instead of step()
        """
        #start = time.time()
        super().propogate(B_override=B_override, step_override=self.__set_stepsize())        
        self.pdf = PDFmodel.pdf(self.position, self.A, self.E, algorithm='simps')
        #print(time.time() - start)
        self.__add_telemetry()


class Incoming(Outgoing):
    pass