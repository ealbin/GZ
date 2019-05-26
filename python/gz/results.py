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

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

from . import coordinates
from . import units

class Result:
    
    def __init__(self, filename, full_telemetry=False):
        self.dirname = os.path.dirname(filename)
        self.filename = os.path.basename(filename)
        self.full_telemetry = full_telemetry
        
        tokens = self.filename.split('_')
        if (not tokens[0].isdecimal()):
            self.theta_mid = int(tokens[1])
            self.phi_mid = int(tokens[2])
        else:
            self.theta_mid = int(tokens[0])
            self.phi_mid = int(tokens[1])
            
    def setZ(self, Z):
        self.Z = Z
    
    def setA(self, A):
        self.A = A
        
    def setE(self, E):
        self.E = E
    
    def setOrigin(self, origin):
        self.origin = origin
        
    def setDist(self, dist):
        self.dist = dist
    
    def setAlgorithm(self, algorithm):
        self.algorithm = algorithm
    
    def setMaxStep(self, max_step):
        self.max_step = max_step
    
    def setRlimit(self, R_limit):
        self.R_limit = R_limit
        
    def setBOverride(self, B_override):
        self.B_override = B_override
        
    def setStepOverride(self, step_override):
        self.step_override = step_override
        
    def setInTelemetry(self, telemetry):
        if (self.full_telemetry):
            self.in_telemetry = telemetry
        else:
            self.in_telemetry = (telemetry[0], telemetry[-1])
        
    def setProtonTelemetry(self, telemetry):
        if (self.full_telemetry):
            self.p_telemetry = telemetry
        else:
            self.p_telemetry = (telemetry[0], telemetry[-1])
        
    def setPDaughterTelemetry(self, telemetry):
        if (self.full_telemetry):
            self.dp_telemetry = telemetry
        else:
            self.dp_telemetry = (telemetry[0], telemetry[-1])
        
    def setNeutronTelemetry(self, telemetry):
        if (self.full_telemetry):
            self.n_telemetry = telemetry
        else:
            self.n_telemetry = (telemetry[0], telemetry[-1])
        
    def setNDaughterTelemetry(self, telemetry):
        if (self.full_telemetry):
            self.dn_telemetry = telemetry
        else:
            self.dn_telemetry = (telemetry[0], telemetry[-1])
        
    def getEarthRadii(self, telemetry):
        pos = telemetry[:3]
        from_earth = pos - coordinates.Cartesian.earth
        dist_Re = np.sqrt(np.dot(from_earth, from_earth)) / (units.SI.radius_earth * units.Change.meter_to_AU)
        return dist_Re
    
    def getSummary(self):
        p_dist  = self.getEarthRadii(self.p_last)
        dp_dist = self.getEarthRadii(self.dp_last)
        n_dist  = self.getEarthRadii(self.n_last)
        dn_dist = self.getEarthRadii(self.dn_last)
        return (p_dist, dp_dist, n_dist, dn_dist)

    def fix(self, telemetry):
        Re = 1.
        
        pos = telemetry[:3] - coordinates.Cartesian.earth
        ro = np.sqrt(np.dot(pos, pos)) / (units.SI.radius_earth * units.Change.meter_to_AU)
        if (ro > 1.):
            sign = -1.
        else:
            sign = 1.
        pos = ro * pos / np.sqrt(np.dot(pos, pos))
        
        beta = telemetry[3:6]
        beta = beta / np.sqrt(np.dot(beta, beta))
        
        ratical = Re*Re - ro*ro + np.dot(pos, beta)**2
        if (ratical < 0.):
            return telemetry[:3]
        
        d = sign * np.dot(pos, beta) + np.sqrt(ratical)
        fixed = pos - (beta * d) 
        fixed = fixed * units.SI.radius_earth * units.Change.meter_to_AU
        fixed = fixed + coordinates.Cartesian.earth
        return fixed
            
    def findLastPos(self):
        p_dist  = self.getEarthRadii(self.p_telemetry[-1])
        dp_dist = self.getEarthRadii(self.dp_telemetry[-1])
        n_dist  = self.getEarthRadii(self.n_telemetry[-1])
        dn_dist = self.getEarthRadii(self.dn_telemetry[-1])
        
        if (p_dist < 1.):
            self.p_last = self.fix(self.p_telemetry[-1])
        else:
            self.p_last = self.p_telemetry[-1][:3]
            
        if (dp_dist < 1.):
            self.dp_last = self.fix(self.dp_telemetry[-1])
        else:
            self.dp_last = self.dp_telemetry[-1][:3]
            
        if (n_dist < 1.):
            self.n_last = self.fix(self.n_telemetry[-1])
        else:
            self.n_last = self.n_telemetry[-1][:3]
            
        if (dn_dist < 1.):
            self.dn_last = self.fix(self.dn_telemetry[-1])
        else:
            self.dn_last = self.dn_telemetry[-1][:3]
            

class Results:
    
    def __init__(self, directory=None, filelist=None, full_telemetry=False):

        if (directory is not None):
            filelist = []
            for file in os.listdir(directory):
                if (file.endswith('.incoming')):
                    filelist.append(os.path.join(directory, file))

        self.directory = directory
        self.filelist = filelist 
        self.results = []

        for file in filelist:
            with open(file, 'r') as f:
                Z = None
                A = None
                E = None
                origin = None
                dist = None
                algorithm = None
                max_step = None
                R_limit = None
                B_override = None
                step_override = None
                in_telemetry = []
                p_telemetry  = []
                dp_telemetry = []
                n_telemetry  = []
                dn_telemetry = []
                seek = 0
                lines = f.readlines()
                for _, line in enumerate(lines):
                    
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
                    
                    # TODO UPDATE TO NEW FORMAT
                    search = '# Origin=['
                    if (line.startswith(search)):
                        tokens = line[len(search):].split()
                        x = float(tokens[0])
                        y = float(tokens[1])
                        z = float(tokens[2].strip(']'))
                        origin = np.asarray([x, y ,z])

                    search = '# Decay_Dist='
                    if (line.startswith(search)):
                        dist = float( line[len(search):].split()[0] )
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
            
                    search = '# Start Telemetry'
                    if (line.startswith(search)):
                        seek = _
                        break
                
                lines = lines[seek + 1:]
                seek = 0
                for _, line in enumerate(lines):
                    search = '#'
                    if (line.startswith(search)):
                        seek = _
                        break
                    in_telemetry.append(np.asarray(line.split(), dtype=np.float64))
                
                if (not lines[seek + 1].startswith('# Proton Telemetry')):
                    print('FORMAT MIS-MATCH')
                    return
                
                lines = lines[seek + 2:]
                seek = 0
                for _, line in enumerate(lines):                 
                    search = '#'
                    if (line.startswith(search)):
                        seek = _
                        break
                    p_telemetry.append(np.asarray(line.split(), dtype=np.float64))
                    
                if (not lines[seek + 1].startswith('# Z-1')):
                    print('FORMAT MIS-MATCH')
                    return
                
                lines = lines[seek + 2:]
                seek = 0
                for _, line in enumerate(lines):                    
                    search = '#'
                    if (line.startswith(search)):
                        seek = _
                        break
                    dp_telemetry.append(np.asarray(line.split(), dtype=np.float64))
                    
                if (not lines[seek + 1].startswith('# Neutron Telemetry')):
                    print('FORMAT MIS-MATCH')
                    return
                
                lines = lines[seek + 2:]
                seek = 0
                for _, line in enumerate(lines):                    
                    search = '#'
                    if (line.startswith(search)):
                        seek = _
                        break
                    n_telemetry.append(np.asarray(line.split(), dtype=np.float64))
                    
                if (not lines[seek + 1].startswith('# Z Daughter')):
                    print('FORMAT MIS-MATCH')
                    return
                
                lines = lines[seek +2:]
                for line in lines:
                    dn_telemetry.append(np.asarray(line.split(), dtype=np.float64))
                    
                result = Result(file, full_telemetry=full_telemetry)
                result.setZ(Z)
                result.setA(A)
                result.setE(E)
                result.setOrigin(origin)
                result.setDist(dist)
                result.setAlgorithm(algorithm)
                result.setMaxStep(max_step)
                result.setRlimit(R_limit)
                result.setBOverride(B_override)
                result.setStepOverride(step_override)
                result.setInTelemetry(in_telemetry)
                result.setProtonTelemetry(p_telemetry)
                result.setPDaughterTelemetry(dp_telemetry)
                result.setNeutronTelemetry(n_telemetry)
                result.setNDaughterTelemetry(dn_telemetry)
                self.results.append(result)

    def HaversineSeparation(pos1, pos2):
        """ Normalized the earth radius aka 1 = Re
        """
        pos1 = np.asarray(pos1)
        pos1 = pos1 / np.sqrt(np.dot(pos1, pos1))
        
        pos2 = np.asarray(pos2)
        pos2 = pos2 / np.sqrt(np.dot(pos2, pos2))
        
        theta1 = np.arccos(pos1[2])
        theta2 = np.arccos(pos2[2])
        
        phi1 = np.arctan2(pos1[1], pos1[0])
        phi2 = np.arctan2(pos2[1], pos2[0])
        
        lat1 = np.pi/2. - theta1
        lat2 = np.pi/2. - theta2
        
        lon1 = phi1
        lon2 = phi2
        
        part1 = np.sin( (lat2 - lat1) / 2. )**2.
        part2 = np.cos(lat1) * np.cos(lat2)
        part3 = np.sin( (lon2 - lon1) / 2. )**2.
    
        return 2. * np.arcsin( np.sqrt(part1 + part2 * part3) )

    def summerize(self, atol=1e-3):
        p_list  = []
        dp_list = []
        n_list  = []
        dn_list = []
        for r in self.results:
            r.findLastPos()
            p, dp, n, dn = r.getSummary()
            p_list.append(p)
            dp_list.append(dp)
            n_list.append(n)
            dn_list.append(dn)
                
        p_near  = np.isclose(p_list,  [1.], atol=atol)
        dp_near = np.isclose(dp_list, [1.], atol=atol)
        n_near  = np.isclose(n_list,  [1.], atol=atol)
        dn_near = np.isclose(dn_list, [1.], atol=atol)
       
        p_dp_both = p_near * dp_near
        n_dn_both = n_near * dn_near
    
        p_xor_dp = (p_near * ~dp_near) + (~p_near * dp_near)
        n_xor_dn = (n_near * ~dn_near) + (~n_near * dn_near)

        p_both = len(p_dp_both[p_dp_both])
        n_both = len(n_dn_both[n_dn_both])
        p_solo = len(p_xor_dp[p_xor_dp])
        n_solo = len(n_xor_dn[n_xor_dn])
        print('N sims: ' + str(len(self.results)))
        print('Proton both:    ' + str(p_both) + ', ' + str(p_both / len(self.results) * 100.) + '%')
        print('One, not both:  ' + str(p_solo) + ', ' + str(p_solo / len(self.results) * 100.) + '%')
        print('Neutron both:   ' + str(n_both) + ', ' + str(n_both / len(self.results) * 100.) + '%')
        print('One, not both:  ' + str(n_solo) + ', ' + str(n_solo / len(self.results) * 100.) + '%')
                
        fig1 = plt.figure(figsize=[15,15])
        bins = np.linspace(0., 1 + 100 * atol, 100)
        plt.hist(p_list,  bins=bins, density=True, log=True, color=mpl.colors.to_rgba('b',.3), label='Proton')
        plt.hist(dp_list, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('m',.3), label='Z-1 Daughter')
        plt.hist(n_list,  bins=bins, density=True, log=True, color=mpl.colors.to_rgba('r',.3), label='Neutron')
        plt.hist(dn_list, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('y',.3), label='Z Daughter')
        plt.legend()
        plt.show()
        
        p_dp_dist = []
        for pair in np.asarray(self.results)[p_dp_both]:
            p_dp_dist.append(Results.HaversineSeparation(pair.p_last, pair.dp_last))
            
        n_dn_dist = []
        for pair in np.asarray(self.results)[n_dn_both]:
            n_dn_dist.append(Results.HaversineSeparation(pair.n_last, pair.dn_last))
        
        p_dp_dist = np.asarray(p_dp_dist) * units.SI.radius_earth
        n_dn_dist = np.asarray(n_dn_dist) * units.SI.radius_earth
        
        fig2 = plt.figure(figsize=[15,15])
        bins = np.logspace(-10,10,1000)
        plt.hist(p_dp_dist, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('b',.3), label='proton-daughter')
        plt.hist(n_dn_dist, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('r',.3), label='neutron-daughter')
        plt.xscale('log')
        plt.legend()
        plt.show()