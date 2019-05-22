#!/usr/bin/env python3

"""Earth
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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from . import coordinates
from . import units

class Patch:
    
    def __init__(self, phi_lo, phi_hi, theta_lo, theta_hi):
        """phi = azimuthal angle = 0 to 360 deg from x axis
        theta = polar angle = 0 to 180 deg from z axis
        """
        self.phi_lo = phi_lo
        self.phi_hi = phi_hi
        self.theta_lo = theta_lo
        self.theta_hi = theta_hi
        
        self.phi_mid = (phi_lo + phi_hi) / 2.
        self.theta_mid = (theta_lo + theta_hi) / 2.
        
        p = self.phi_mid * (np.pi / 180.)
        t = self.theta_mid * (np.pi / 180.)
        x = np.sin(t) * np.cos(p)
        y = np.sin(t) * np.sin(p)
        z = np.cos(t)
        self.zenith = np.asarray([x, y, z])
        
        
class Earth:
    
    OUT_JOB_PATH = './out_jobs'
    
    def run(wedges=4, bands=3, Zlist=[2, 26, 92], Elist=[2e18, 20e18, 200e18], runs=1000):
        earth = Earth()
        for z in Zlist:
            for e in Elist:
                earth.outgoing_jobs(z, e, max_step=.01, runs=runs, cone=90.)
    
    def __init__(self, wedges=4, bands=3):
        self.wedges = wedges
        self.bands = bands
        
        self.phi_sep = 360. / wedges
        self.theta_sep = 180. / bands
        
        self.phi_offset = self.phi_sep / 2.
        self.theta_offset = 0.        
        
        self.patches = []
        for w in range(wedges):
            for b in range(bands):
                phi_lo = self.phi_offset + w * self.phi_sep
                phi_hi = phi_lo + self.phi_sep
                theta_lo = self.theta_offset + b * self.theta_sep
                theta_hi = theta_lo + self.theta_sep
                self.patches.append(Patch(phi_lo, phi_hi, theta_lo, theta_hi))


    def draw(self, ax=None):
        if (ax is None):
            fig = plt.figure(figsize=[16,16])
            ax = plt.axes(projection='3d')
        
        for patch in self.patches:      
            phi_lo = patch.phi_lo * np.pi / 180.
            phi_hi = patch.phi_hi * np.pi / 180.
            theta_lo = patch.theta_lo * np.pi / 180.
            theta_hi = patch.theta_hi * np.pi / 180.
            
            u, v = np.mgrid[phi_lo:phi_hi:10j, theta_lo:theta_hi:10j]
            r = units.SI.radius_earth * units.Change.meter_to_AU
            x = r * np.cos(u)*np.sin(v)
            y = r * np.sin(u)*np.sin(v)
            z = r * np.cos(v)
            x += coordinates.Cartesian.earth[0]
            y += coordinates.Cartesian.earth[1]
            z += coordinates.Cartesian.earth[2]
            ax.plot_surface(x, y, z, color=tuple(np.random.rand(3)))

 
    def outgoing_jobs(self, Z, E, max_step=.01, A=None, R_limit=None, runs=100, cone=90., 
                      seed=None, out_path=None, job_path=None, name_header=None, name_tail=None,
                      B_override=None):
        if (seed is not None):
            np.random.seed(seed)
        
        if (job_path is None):
            job_path = Earth.OUT_JOB_PATH
        
        if (not os.path.isdir(job_path)):
            os.makedirs(job_path)
        
        for patch in self.patches:
            p_mid = int(patch.phi_mid)
            t_mid = int(patch.theta_mid)

            if (name_header is not None):
                filename = name_header + '_'
            else:
                filename = ''
            filename += str(t_mid) + '_' + str(p_mid)
            
            if (name_tail is not None):
                filename += '_' + name_tail
            else:
                filename += '_' + str(Z) + '_' + str(int(E/1e18))    
            filename += '.py'

            position = coordinates.Cartesian.earth
            position = position + patch.zenith * units.SI.radius_earth * units.Change.meter_to_AU

            with open(os.path.join(job_path, filename), 'w') as f:
                f.write('#!/usr/bin/env python3\n')
                f.write('#\n')
                f.write('# Outgoing propagation job: ' + __version__ + '\n')
                f.write('# Time of writing: ' + str(datetime.datetime.now()) + '\n')        
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
                f.write('# Setup\n')
                f.write('# Wedges=' + str(self.wedges) + '\n')
                f.write('# Bands=' + str(self.bands) + '\n')
                f.write('# Runs=' + str(runs) + '\n')
                f.write('# Cone=' + str(cone) + ' [deg]\n')
                f.write('# Seed=' + str(seed) + '\n')
                f.write('#\n')
                f.write('# Patch\n')
                f.write('# Phi_lo=' + str(patch.phi_lo) + ' [deg]\n')
                f.write('# Phi_mid=' + str(patch.phi_mid) + ' [deg]\n')
                f.write('# Phi_hi=' + str(patch.phi_hi) + ' [deg]\n')
                f.write('# Theta_lo=' + str(patch.theta_lo) + ' [deg]\n')
                f.write('# Theta_mid=' + str(patch.theta_mid) + ' [deg]\n')
                f.write('# Theta_hi=' + str(patch.theta_hi) + ' [deg]\n')
                f.write('# Zenith=' + str(patch.zenith) + '\n')
                f.write('#\n')
                f.write('# Parameters\n')
                f.write('# Z=' + str(np.abs(Z)) + ' [proton number]\n')
                f.write('# A=' + str(A) + ' [atomic mass units]\n')
                f.write('# E=' + str(E) + ' [electron volts]\n')
                f.write('# Max_step=' + str(max_step) + ' [AU]\n')
                f.write('# R_limit=' + str(R_limit) + ' [AU]\n')
                f.write('#\n')
                f.write('# Script\n\n')
                f.write('import gz\n\n')
                
                phis = 2.*np.pi * np.random.random(runs)
                thetas = cone * np.pi / 180. * np.random.random(runs)
                for t, p in zip(thetas, phis):
                    r  = np.cos(t)
                    th = np.sin(t) * np.cos(p)
                    ph = np.sin(t) * np.sin(p)
                    theta = np.arccos(patch.zenith[2])
                    phi   = np.arctan2(patch.zenith[1], patch.zenith[0])
                    beta = coordinates.Spherical.toCartesian([r, th, ph], theta, phi)
                    
                    if (name_header is not None):
                        out_name = name_header + '_'
                    else:
                        out_name = ''
                    out_name += str(t_mid) + '_' + str(p_mid) + '_'
                    out_name += str(Z) + '_' + str(int(E/1e18)) + '_'
                    out_name += str(int(t * 180. / np.pi)) + '_' + str(int(p * 180./np.pi))
                    
                    args  = '[' + str(position[0]) + ', ' + str(position[1]) + ', ' + str(position[2]) + '], '
                    args += '[' + str(beta[0])     + ', ' + str(beta[1])     + ', ' + str(beta[2])     + '], '
                    args += str(Z) + ', '
                    args += str(E) + ', '
                    args += 'A=' + str(A) + ', '
                    args += 'max_step=' + str(max_step) + ', '
                    if (R_limit is not None):
                        args += 'R_limit=' + str(R_limit) + ', '
                    args += 'save_path=' + str(out_path) + ', '
                    args += 'filename=' + "'" + out_name + "'"
                    f.write('outgoing = gz.path.Outgoing(' + args + ')\n') 
                    args = ''
                    if (B_override is not None):
                        b_str = '[' + str(B_override[0]) + ', ' + str(B_override[1]) + ', ' + str(B_override[2]) + ']'
                        args = 'B_override=' + b_str                    
                    f.write('outgoing.propagate(' + args + ')\n\n')
      
        
    def incoming_jobs(directory=None, filelist=None):

        if (directory is not None):
            filelist = []
            for file in os.listdir(directory):
                if (file.endswith('.outgoing')):
                    filelist.append(file)

        for file in filelist:
            with open(file, 'r') as f:
                Z = None
                A = None
                E = None
                algorithm = None
                max_step = None
                R_limit = None
                B_override = None
                step_override = None
                telemetry = []
                seek = 0
                for _, line in enumerate(f.readlines()):
                    
                    search = '# Z='
                    if (line.startswith(search)):
                        Z = int( line[len(search):].split()[0] )
                        continue
                    
                    search = '# A='
                    if (line.startswith(search)):
                        A = line[len(search):].split()[0]
                        try:
                            A = float(A)
                        except ValueError:
                            A = None
                        continue
                    
                    search = '# E='
                    if (line.startswith(search)):
                        E = float( line[len(search):].split()[0] )
                        continue
                    
                    search = '# Algorithm='
                    if (line.startswith(search)):
                        algorithm = line[len(search):].split()[0]
                        continue                    

                    search = '# Max_Step='
                    if (line.startswith(search)):
                        max_step = line[len(search):].split()[0]
                        try:
                            max_step = float(max_step)
                        except ValueError:
                            max_step = None
                        continue

                    search = '# R_Limit='
                    if (line.startswith(search)):
                        R_limit = line[len(search):].split()[0]
                        try:
                            R_limit = float(R_limit)
                        except ValueError:
                            R_limit = None
                        continue
            
                    search = '# B_Override='
                    if (line.startswith(search)):
                        B_override = line[len(search):].split()
                        try:
                            B_override = np.asarray(B_override[:3], dtype=np.float64)
                        except ValueError:
                            B_override = None
                        continue                        

                    search = '# Step_Override='
                    if (line.startswith(search)):
                        step_override = line[len(search):].split()[0]
                        try:
                            step_override = float(step_override)
                        except ValueError:
                            step_override = None
                        continue
            
                    search = '# Telemetry'
                    if (line.startswith(search)):
                        f.seek(0)
                        seek = _
                        break
                
                for line in f.readlines()[seek + 1:]:
                    telemetry.append(np.asarray(line.split(), dtype=np.float64))
                
                origin = telemetry[0][:3]
                position = telemetry[-1][:3]
                beta = -1 * telemetry[-1][3:6]
                
                print('max step = ' + str(max_step))


#    def __init__(self, origin, position, beta, Z, A, E, decay_dist, 
#                 max_step=None, R_limit=6., save=True, save_path=None, filename=None):