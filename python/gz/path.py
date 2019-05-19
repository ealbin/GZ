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

import datetime
import numpy as np
import os
import platform
import time

from scipy import integrate

from . import coordinates
from . import magnetic_field
from . import units

class Path:
    
    EULER_DIVISOR  = 1e4
    DOP853_DIVISOR = 1e2
    
    def __init__(self, position, beta, Z, E, max_step=1.):
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
        
        self.max_step = max_step
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
            self.step = self.max_step
        else:
            B = np.sqrt(np.dot(self.B, self.B))
            gyro_radius  = 1. / units.SI.lightspeed / np.abs(self.ratio) / B
            gyro_radius *= units.Change.meter_to_AU
            self.step    = min(self.max_step, gyro_radius / self.step_divisor)
            
    
    def propagate(self, B_override=None, step_override=None, algorithm='dop853'):
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
    
    LIMIT_BUFFER  = 2. # AU
    
    DEFAULT_SAVE_PATH = './telemetry'
    
    def __init__(self, position, beta, Z, E,
                 A=None, max_step=None, R_limit=6., zigzag=False, 
                 save=True, save_path=None, filename=None):
        """
        position: np.array(x,y,z) AU start position on earth
        beta: np.array(bx,by,bz) start beta (away from earth)
        Z: atomic number
        A: atomic mass if none then auto assign
        E: energy in eV
        R_limit: radius [AU] of maximum propagation
        """        
        if (max_step is None):
            Path.__init__(self, position, beta, -Z, E)
        else:
            Path.__init__(self, position, beta, -Z, E, max_step=max_step)
        
        self.A = A
        self.R_limit = R_limit
        
        self.telementry = [np.concatenate([self.position, self.beta, [self.distance]])]
        self.last_save = self.distance
        self.save_distance = self.max_step / 10.
        
        self.zigzag = zigzag
        self.save = save
        self.save_path = save_path
        self.filename = filename
         
    def _add_telemetry(self):
        near_sun   = units.SI.radius_sun   * 10. * units.Change.meter_to_AU
        near_earth = units.SI.radius_earth * 2.  * units.Change.meter_to_AU
        
        if (self.dist_sun < near_sun or self.dist_earth < near_earth
            or self.distance - self.last_save >= self.save_distance):
            
            self.telementry.append(np.concatenate([self.position, self.beta, [self.distance]]))     
            self.last_save = self.distance
           
    def _set_stepsize(self):
        # any special needs here
        return Path._set_stepsize(self)
    
    def propagate(self, B_override=None, step_override=None, algorithm='dop853'):
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
        
        start = time.time()
        while (stop_condition()):
            Path.propagate(self, B_override=self.B, step_override=self.step, algorithm=algorithm)        
            self._add_telemetry()
        self.elapsed_sec = time.time() - start
        
        if (self.last_save < self.distance):
            self.telementry.append(np.concatenate([self.position, self.beta, [self.distance]]))
        
        if (self.save):
            self.algorithm = algorithm
            self.save_telemetry()
        
    def save_telemetry(self):
        if (self.save_path is None):
            self.save_path = Outgoing.DEFAULT_SAVE_PATH
        
        if (not os.path.isdir(self.save_path)):
            os.makedirs(self.save_path)
        
        if (self.filename is None):
            self.filename  = str(np.abs(self.Z))
            self.filename += '_'
            self.filename += str(int(self.E / 1e18))
    
        test_name = self.filename
        full_path = os.path.join(self.save_path, test_name + '.outgoing')
        _ = 1
        while (os.path.exists(full_path)):
            test_name = self.filename + '_' + str(_)
            full_path = os.path.join(self.save_path, test_name + '.outgoing')
            _ += 1
        self.filename = test_name + '.outgoing'
        
        with open(os.path.join(self.save_path, self.filename), 'w') as f:
            f.write('# Outgoing propagation: ' + __version__ + '\n')
            f.write('# Time of writing: ' + str(datetime.datetime.now()) + '\n')        
            f.write('# Run time [sec]: ' + str(self.elapsed_sec) + '\n')
            f.write('#\n')
            f.write('# Platform\n')
            uname = platform.uname()
            f.write('# Node=' + uname.node + '\n')
            f.write('# Machine=' + uname.machine + '\n')
            f.write('# System=' + uname.system + '\n')
            f.write('# Version=' + uname.version + '\n')
            f.write('# Release=' + uname.release + '\n')
            f.write('# Processor=' + uname.processor + '\n')
            f.write('#\n')
            f.write('# Parameters\n')
            f.write('# Z=' + str(np.abs(self.Z)) + ' [proton number]\n')
            f.write('# A=' + str(self.A) + ' [atomic mass units]\n')
            f.write('# E=' + str(self.E) + ' [electron volts]\n')
            f.write('# Algorithm=' + self.algorithm + '\n')
            f.write('#\n')
            f.write('# Key\n')
            f.write('# position_x, position_y, position_z, beta_x, beta_y, beta_z, path_distance\n')
            f.write('# units: positions=AU, beta=unitless, distance=AU\n')
            f.write('#\n')
            f.write('# Telemetry\n')
            for _ in self.telementry:
                for val in _:
                    f.write(str(val) + ' ')
                f.write('\n')

        # pdf stuff going to post-processing
        #self.pdf = PDFmodel.pdf(self.position, self.A, self.E, algorithm='simps')

class Incoming(Outgoing):
    pass