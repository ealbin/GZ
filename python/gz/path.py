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
from . import probability
from . import relativity
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
        if (self.ratio == 0. or np.sqrt(np.dot(self.B, self.B)) == 0.):
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
        
        self.telemetry = [np.concatenate([self.position, self.beta, [self.distance]])]
        self.last_save = self.distance
        self.save_distance = self.max_step / 10.
        
        self.zigzag = zigzag
        self.save = save
        self.save_path = save_path
        self.filename = filename
         
    def _add_telemetry(self):
        near_sun   = units.SI.radius_sun   * 10. * units.Change.meter_to_AU
        near_earth = units.SI.radius_earth * 10. * units.Change.meter_to_AU
        
        if (self.dist_sun < near_sun or self.dist_earth < near_earth
            or self.distance - self.last_save >= self.save_distance):
            
            self.telemetry.append(np.concatenate([self.position, self.beta, [self.distance]]))     
            self.last_save = self.distance
           
    def _set_B(self, B_override=None):
        if (B_override is not None):
            self.B = np.asarray(B_override, dtype=np.float64)
        else:
            if (self.interpolate_B):
                self.B = magnetic_field.cartesianTesla(self.position)         
            else:
                self.B = magnetic_field.cartesianTesla(self.position, close2sun=100.)
        
    def _set_stepsize(self):
        # any special needs here
        Path._set_stepsize(self)
    
    def propagate(self, B_override=None, interpolate_B=True, step_override=None, algorithm='dop853'):
        """
        Propagates one step
        B_override: use this B instead of Bfield
        step_override: use this step instead of step()
        """
        self.B_override = B_override
        self.interpolate_B = interpolate_B
        self.step_override = step_override
        
        self._set_B(B_override=B_override)
        
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
                if (self.distance + self.step > self.R_limit + Outgoing.LIMIT_BUFFER):
                    self.step = self.R_limit + Outgoing.LIMIT_BUFFER - self.distance
                return self.distance < self.R_limit + Outgoing.LIMIT_BUFFER
        else:
            def stop_condition():
                return self.distance < self.R_limit + Outgoing.LIMIT_BUFFER and self.dist_sun < self.R_limit
        
        start = time.time()
        while (stop_condition()):
            Path.propagate(self, B_override=self.B, step_override=self.step, algorithm=algorithm) 
            self._set_B(B_override=B_override)
            self._add_telemetry()
            if (step_override is None):
                self._set_stepsize()
        self.elapsed_sec = time.time() - start
        
        if (self.last_save < self.distance):
            self.telemetry.append(np.concatenate([self.position, self.beta, [self.distance]]))
        
        if (self.save):
            self.algorithm = algorithm
            self.save_telemetry()
  
    def save_telemetry(self):
        if (self.save_path is None):
           self.save_path = Outgoing.DEFAULT_SAVE_PATH
            
        subdir = str(np.abs(self.Z)) + '_' + str(int(self.E/1e15))
        self.save_path = os.path.join(self.save_path, subdir)   
        if (not os.path.isdir(self.save_path)):
            os.makedirs(self.save_path)
        
        if (self.filename is None):
            self.filename  = str(np.abs(self.Z))
            self.filename += '_'
            self.filename += str(int(self.E / 1e15))
    
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
            f.write('# Max_Step=' + str(self.max_step) + ' [AU]\n')
            f.write('# R_Limit=' + str(self.R_limit) + ' [AU]\n')
            B_str = str(self.B_override)
            if (self.B_override is not None):
                B_str = str(self.B_override[0]) + ' ' + str(self.B_override[1]) + ' ' + str(self.B_override[2])
            f.write('# B_Override=' + B_str + ' [T]\n')
            f.write('# Step_Override=' + str(self.step_override) + '\n')
            f.write('#\n')
            f.write('# Key\n')
            f.write('# position_x, position_y, position_z, beta_x, beta_y, beta_z, path_distance\n')
            f.write('# units: positions=AU, beta=unitless, distance=AU\n')
            f.write('#\n')
            f.write('# Telemetry\n')
            for _ in self.telemetry:
                for val in _:
                    f.write(str(val) + ' ')
                f.write('\n')


class Incoming(Outgoing):
    
    def __init__(self, origin, position, beta, Z, A, E, decay_dist, 
                 max_step=None, R_limit=6., save=True, save_path=None, filename=None):
        
        if (max_step is None):
            Path.__init__(self, position, beta, Z, E)
        else:
            Path.__init__(self, position, beta, Z, E, max_step=max_step)
        
        self.origin = np.asarray(origin)
        self.decay_dist = decay_dist
        self.R_limit = R_limit
        
        if (A is None):
            self.A = units.Nuclide.mass_number(Z)
        else:
            self.A = A
        
        self.telemetry = [np.concatenate([self.position, self.beta, [self.distance]])]
        self.last_save = self.distance
        self.save_distance = self.max_step / 10.
        self.near_earth = False
        
        self.save = save
        self.save_path = save_path
        self.filename = filename
        
    def _add_telemetry(self):
        # add anything custom
        Outgoing._add_telemetry(self)
            
    def _set_stepsize(self):
        if (self.near_earth or self.dist_earth <= self.max_step):
            self.near_earth = True
            if (self.dist_earth > self.max_step):
                self.near_earth = False
            self.save_distance = 10. * units.SI.radius_earth * units.Change.meter_to_AU
            if (self.dist_earth > 2 * units.SI.radius_earth * units.Change.meter_to_AU):
                self.step = units.SI.radius_earth * units.Change.meter_to_AU / 5.
            else:
                self.step = units.SI.radius_earth * units.Change.meter_to_AU / 50.
        else:
            Path._set_stepsize(self)         
            
    def propagate(self, B_override=None, step_override=None, algorithm='dop853', seed=None):
        """
        Propagates one step
        B_override: use this B instead of Bfield
        step_override: use this step instead of step()
        """
        self.B_override = B_override
        self.step_override = step_override

        Outgoing._set_B(self, B_override=B_override)
        
        if (step_override is not None):
            self.step = step_override
        else:
            if (algorithm == 'euler'):
                self.step_divisor = Path.EULER_DIVISOR
            elif (algorithm == 'dop853'):
                self.step_divisor = Path.DOP853_DIVISOR
            Path._set_stepsize(self)
            self._set_stepsize()        
        
        def remaining():
            return self.decay_dist - self.distance
     
        # Propogate nucleus until time to dissintegrate
        start = time.time()
        while (remaining() > 0):
            if (remaining() < self.step):
                self.step = remaining()
            Path.propagate(self, B_override=self.B, step_override=self.step, algorithm=algorithm)  
            Outgoing._set_B(self, B_override=B_override)
            self._add_telemetry()
            if (step_override is None):
                self._set_stepsize()
        if (self.last_save < self.distance):
            self.telemetry.append(np.concatenate([self.position, self.beta, [self.distance]]))

        # Photodissintegration
        # "1" = original nucleus
        # "2" = solar photon
        # "3" = proton or neutron
        # "4" = daughter nucleus
        e1 = self.E
        p1 = relativity.momentum(e1, self.A * units.Change.amu_to_eV, self.beta)
        e2 = probability.Solar.get_photon(self.position, self.beta, self.Z, self.E, seed=seed) # seed is set here
        p2 = relativity.momentum(e2, 0., self.position / self.dist_sun)
        
        Epn =  1.          / self.A * (self.E + e2) # proton/neutron energy
        Ed  = (self.A - 1) / self.A * (self.E + e2) # daugter nucleous energy
        Zp  = 1 # proton charge
        Zn  = 0 # neutron charge
        Zdp = self.Z - 1 # daughter (proton ejection) charge
        Zdn = self.Z     # daughter (neutron ejection) charge

        e3 = Epn
        m3 = 1. * units.Change.amu_to_eV
        e4 = Ed
        m4 = (self.A - 1) * units.Change.amu_to_eV
        
        # "p" is the net 3-momentum
        p       = p1 + p2
        p_mag   = np.sqrt(np.dot(p, p))
        p_hat   = p / p_mag
        p_theta = np.arccos(p_hat[2])
        p_phi   = np.arctan2(p_hat[1], p_hat[0])
        
        p3_mag = relativity.momentum_mag(e3, m3)
        p4_mag = relativity.momentum_mag(e4, m4)
        
        # Angle between p3 and p4
        theta = relativity.theta(e1, p1, e2, p2, e3, m3, e4, m4)
        # Angle between p3 and p
        cosTheta = (p3_mag * p4_mag * np.cos(theta) + p3_mag**2) / (p3_mag * p_mag)
        if (cosTheta > 1. and np.isclose(cosTheta, 1.)):
            cosTheta = 1.
        if (cosTheta < -1. and np.isclose(cosTheta, -1.)):
            cosTheta = -1.
        theta3 = np.arccos(cosTheta)
        # Azimuthal angle around p
        phi3 = 2. * np.pi * np.random.random()
        
        p3_r = p3_mag * np.cos(theta3)
        p3_t = p3_mag * np.sin(theta3) * np.cos(phi3)
        p3_p = p3_mag * np.sin(theta3) * np.sin(phi3) 
        p3 = coordinates.Spherical.toCartesian(np.asarray([p3_r, p3_t, p3_p]), p_theta, p_phi)
        p4 = p - p3
        
        beta_3 = p3 / p3_mag
        beta_4 = p4 / p4_mag
        
        self.p_path  = Incoming(None, self.position, beta_3, Zp,  None, e3, None, max_step=self.max_step) # ejected proton
        self.n_path  = Incoming(None, self.position, beta_3, Zn,  None, e3, None, max_step=self.max_step) # ejected neutron
        self.dp_path = Incoming(None, self.position, beta_4, Zdp, None, e4, None, max_step=self.max_step) # Z-1 nucleus
        self.dn_path = Incoming(None, self.position, beta_4, Zdn, None, e4, None, max_step=self.max_step) # A-1 nucleus

        for subpath in [self.p_path, self.n_path, self.dp_path, self.dn_path]:
            subpath.sub_propogate(B_override=B_override, step_override=None, algorithm=algorithm)
                        
        self.elapsed_sec = time.time() - start
            
        if (self.save):
            self.algorithm = algorithm
            self.save_telemetry()
            
    # Sub-propogate children
    def sub_propogate(self, B_override=None, step_override=None, algorithm='dop853'):

        self.B_override = B_override
        self.step_override = step_override

        Outgoing._set_B(self, B_override=B_override)
        
        if (step_override is not None):
            self.step = step_override
        else:
            if (algorithm == 'euler'):
                self.step_divisor = Path.EULER_DIVISOR
            elif (algorithm == 'dop853'):
                self.step_divisor = Path.DOP853_DIVISOR
            Path._set_stepsize(self)
            self._set_stepsize()        
        
        dist_earth_init = self.dist_earth
        
        def keep_going():
            if (self.dist_earth > (1.01) * units.SI.radius_earth * units.Change.meter_to_AU
                and self.distance < dist_earth_init + Outgoing.LIMIT_BUFFER):
                return True
            return False

        while (keep_going()):
            Path.propagate(self, B_override=self.B, step_override=self.step, algorithm=algorithm)
            Outgoing._set_B(self, B_override=B_override)
            self._add_telemetry()  
            if (step_override is None):
                self._set_stepsize()            
        if (self.last_save < self.distance):
            self.telemetry.append(np.concatenate([self.position, self.beta, [self.distance]]))
            
            
    def save_telemetry(self):
        if (self.save_path is None):
            self.save_path = Outgoing.DEFAULT_SAVE_PATH
        
        self.save_path = os.path.join(self.save_path, str(self.Z))   
        if (not os.path.isdir(self.save_path)):
            os.makedirs(self.save_path)
        
        if (self.filename is None):
            self.filename  = str(np.abs(self.Z))
            self.filename += '_'
            self.filename += str(int(self.E / 1e15))
    
        test_name = self.filename
        full_path = os.path.join(self.save_path, test_name + '.incoming')
        _ = 1
        while (os.path.exists(full_path)):
            test_name = self.filename + '_' + str(_)
            full_path = os.path.join(self.save_path, test_name + '.incoming')
            _ += 1
        self.filename = test_name + '.incoming'
        
        with open(os.path.join(self.save_path, self.filename), 'w') as f:
            f.write('# Incoming propagation: ' + __version__ + '\n')
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
            f.write('# Origin=' + str(self.origin[0]) + ' ' + str(self.origin[1]) + ' ' + str(self.origin[2]) + ' [AU]\n')
            f.write('# Decay_Dist=' + str(self.decay_dist) + ' [AU]\n')
            f.write('# Algorithm=' + self.algorithm + '\n')
            f.write('# Max_Step=' + str(self.max_step) + ' [AU]\n')
            f.write('# R_Limit=' + str(self.R_limit) + ' [AU]\n')
            B_str = str(self.B_override)
            if (self.B_override is not None):
                B_str = str(self.B_override[0]) + ' ' + str(self.B_override[1]) + ' ' + str(self.B_override[2])
            f.write('# B_Override=' + B_str + ' [T]\n')
            f.write('# Step_Override=' + str(self.step_override) + '\n')
            f.write('#\n')
            f.write('# Key\n')
            f.write('# position_x, position_y, position_z, beta_x, beta_y, beta_z, path_distance\n')
            f.write('# units: positions=AU, beta=unitless, distance=AU\n')
            f.write('#\n')
                    
            f.write('# Start Telemetry\n')
            for _ in self.telemetry:
                for val in _:
                    f.write(str(val) + ' ')
                f.write('\n')
            f.write('#\n')
                    
            f.write('# Proton Telemetry\n')
            for _ in self.p_path.telemetry:
                for val in _:
                    f.write(str(val) + ' ')
                f.write('\n')
            f.write('#\n')
                    
            f.write('# Z-1 Daughter Telemetry\n')
            for _ in self.dp_path.telemetry:
                for val in _:
                    f.write(str(val) + ' ')                    
                f.write('\n')
            f.write('#\n')

            f.write('# Neutron Telemetry\n')
            for _ in self.n_path.telemetry:
                for val in _:
                    f.write(str(val) + ' ')
                f.write('\n')
            f.write('#\n')
                    
            f.write('# Z Daughter Telemetry\n')
            for _ in self.dn_path.telemetry:
                for val in _:
                    f.write(str(val) + ' ')                    
                f.write('\n')