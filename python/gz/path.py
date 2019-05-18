#!/usr/bin/env python3

"""
Description
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

from scipy import integrate

from . import coordinates
from . import magnetic_field
from . import probability
from . import units

class Path:
    
    MAX_STEP = 1. # AU
    EULER_DIVISOR  = 1e4
    DOP853_DIVISOR = 1e2
    
    def __init__(self, position, beta, Z, E):
        """
        position: np.array(x,y,z) start position
        beta: np.array(bx,by,bz) start beta (direction of propagation)
        Z: atomic number
        E: energy in eV
        """
        self.position = np.asarray(position, dtype=np.float64)
        self.beta = np.asarray(beta, dtype=np.float64)
        self.beta = self.beta / np.sqrt(np.dot(self.beta, self.beta))
        self.Z = Z
        self.E = E
        self.ratio = Z / E
        
        self.distance = 0.
        self._set_dist_earth()
        self._set_dist_sun()
        
    def _set_dist_earth(self):
        r = self.position - coordinates.Cartesian.earth
        self.dist_earth = np.sqrt(np.dot(r, r))
        
    def _set_dist_sun(self):
        r = self.position - coordinates.Cartesian.sun
        self.dist_sun = np.sqrt(np.dot(r, r))

    def _set_stepsize(self):
        if (self.ratio == 0.):
            self.step = Path.MAX_STEP
        else:
            B = np.sqrt(np.dot(self.B, self.B))
            gyro_radius  = 1. / units.SI.lightspeed / np.abs(self.ratio) / B
            gyro_radius *= units.Change.meter_to_AU
            self.step    = min(Path.MAX_STEP, gyro_radius / self.step_divisor)
            
    
    def propogate(self, B_override=None, step_override=None, algorithm='dop853'):
        """
        Propagates one step
        B_override: use this B instead of Bfield [tesla]
        step_override: use this step instead of step()
        """        
        if (B_override is not None):
            self.B = np.asarray(B_override, dtype=np.float64)
        else:
            self.B = magnetic_field.cartesianTesla(self.position) 
        
        if (step_override is not None):
            self.step = step_override
        else:
            if (algorithm == 'euler'):
                self.step_divisor = Path.EULER_DIVISOR
            elif (algorithm == 'dop853'):
                self.step_divisor = Path.DOP853_DIVISOR
            self._set_stepsize()        
            
        if (algorithm == 'euler'):
            dbeta_ds  = units.Change.AU_to_meter
            dbeta_ds *= self.ratio
            dbeta_ds *= np.cross(self.beta, units.SI.lightspeed * self.B)
    
            self.beta += dbeta_ds * self.step
            self.beta = self.beta / np.sqrt(np.dot(self.beta, self.beta))    
            self.position += self.beta * self.step

        else:
            def ode(t, Y):
                beta = Y[3:]
                dbeta_ds  = units.Change.AU_to_meter
                dbeta_ds *= self.ratio
                dbeta_ds *= np.cross(beta, units.SI.lightspeed * self.B)
                return np.concatenate([beta, dbeta_ds])
            try:
                self.integrator
            except (AttributeError, NameError):
                if (algorithm == 'dop853'):
                    self.integrator = integrate.ode(ode).set_integrator('dop853')    

            initial_conditions = np.concatenate([self.position, self.beta])
            self.integrator.set_initial_value(initial_conditions, 0.)            
            self.integrator.integrate(self.integrator.t + self.step)
            self.position = self.integrator.y[:3]
            self.beta = self.integrator.y[3:]
            self.beta = self.beta / np.sqrt(np.dot(self.beta, self.beta))

        self.distance += self.step
        self._set_dist_earth()
        self._set_dist_sun()
    
    
class Outgoing(Path):
    
    SAVE_DISTANCE = Path.MAX_STEP / 10. # AU
    LIMIT_BUFFER  = 2. # AU
    
    def __init__(self, position, beta, Z, E, A=None, R_limit=6., zigzag=False):
        """
        position: np.array(x,y,z) AU start position on earth
        beta: np.array(bx,by,bz) start beta (away from earth)
        Z: atomic number
        A: atomic mass if none then auto assign
        E: energy in eV
        R_limit: radius [AU] of maximum propagation
        """        
        super().__init__(position, beta, -Z, E)
        self.A = A
        self.R_limit = R_limit
        
        self.telementry = [np.concatenate([self.position, self.beta, [self.distance]])]
        self.last_save = self.distance
        
        self.zigzag = zigzag
         
    def _add_telemetry(self):
        near_sun = units.SI.radius_sun * 10. * units.Change.meter_to_AU
        near_earth = units.SI.radius_earth * 2. * units.Change.meter_to_AU
        if (self.dist_sun < near_sun or self.dist_earth < near_earth
            or self.distance - self.last_save >= Outgoing.SAVE_DISTANCE):
            self.telementry.append(np.concatenate([self.position, self.beta, [self.distance]]))     
            self.last_save = self.distance
           
    def _set_stepsize(self):
        # any special needs here
        return super()._set_stepsize()
    
    def propogate(self, B_override=None, step_override=None, algorithm='dop853'):
        """
        Propagates one step
        B_override: use this B instead of Bfield
        step_override: use this step instead of step()
        """
        if (B_override is not None):
            self.B = np.asarray(B_override, dtype=np.float64)
        else:
            self.B = magnetic_field.cartesianTesla(self.position) 
        
        if (step_override is not None):
            self.step = step_override
        else:
            if (algorithm == 'euler'):
                self.step_divisor = Path.EULER_DIVISOR
            elif (algorithm == 'dop853'):
                self.step_divisor = Path.DOP853_DIVISOR
            self._set_stepsize()        
        
        if (self.zigzag):
            def stop_condition():
                return self.distance < self.R_limit + Outgoing.LIMIT_BUFFER
        else:
            def stop_condition():
                return self.distance < self.R_limit + Outgoing.LIMIT_BUFFER and self.dist_sun < self.R_limit
        while (stop_condition()):
            super().propogate(B_override=self.B, step_override=self.step, algorithm=algorithm)        
            self._add_telemetry()
            
        if (self.last_save < self.distance):
            self.telementry.append(np.concatenate([self.position, self.beta, [self.distance]]))
        
    def save_telemetry(self):
        pass
        # check if telemetry folder exists
        # make filename
        # write header
        # write telementry

        # pdf stuff going to post-processing
        #self.pdf = PDFmodel.pdf(self.position, self.A, self.E, algorithm='simps')

class Incoming(Outgoing):
    pass