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
        if (tokens[0][0].isalpha()):
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
        Re = units.SI.radius_earth * units.Change.meter_to_AU
        
        pvec = telemetry[:3] - coordinates.Cartesian.earth
        pmag = np.sqrt(np.dot(pvec, pvec))
        phat = pvec / pmag
        if (pmag == Re):
            return telemetry[:3]
        
        bhat = telemetry[3:6]
        bhat = bhat / np.sqrt(np.dot(bhat,bhat))
        
        cosTheta = np.dot(phat, bhat)
        discriminant = Re*Re - pmag*pmag * (1. - cosTheta*cosTheta)    
        if (discriminant < 0.):
            return telemetry[:3]
        
        cplus  = -1. * pmag * cosTheta + np.sqrt(discriminant)
        cminus = -1. * pmag * cosTheta - np.sqrt(discriminant)
        c = np.asarray([cplus, cminus])
        
        if (pmag > Re):
            if (cosTheta > 0.):
                c = c[c < 0.]
                c = -1. * np.min(np.abs(c))                
            else:
                c = c[c > 0.]
                c = np.min(c)
        else:
            c = c[c < 0.]
            if (cosTheta > 0.):
                c = -1. * np.max(np.abs(c))
            else:
                c = -1. * np.min(np.abs(c))
        
        pvec = (pvec + c * bhat) + coordinates.Cartesian.earth
        return pvec
            
    def findLastPos(self):
        self.p_last  = self.fix(self.p_telemetry[-1])
        self.dp_last = self.fix(self.dp_telemetry[-1])
        self.n_last  = self.fix(self.n_telemetry[-1])
        self.dn_last = self.fix(self.dn_telemetry[-1]) 
        return
        
            

class Results:
    
    def __init__(self, directory=None, filelist=None, full_telemetry=False, cone=None):

        if (directory is not None):
            filelist = []
            for file in os.listdir(directory):
                if (file.endswith('.incoming')):
                    if (cone is not None):
                        tokens = file.split('_')
                        if (tokens[0][0].isalpha()):
                            tokens = tokens[1:]
                        theta = float(tokens[4])
                        if (theta < cone):
                            filelist.append(os.path.join(directory, file))
                    else:
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

    def summerize(self, atol=1e-2):
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

        p_neither = ~(p_dp_both + p_xor_dp)
        n_neither = ~(n_dn_both + n_xor_dn)
        
        p_both = len(p_dp_both[p_dp_both])
        n_both = len(n_dn_both[n_dn_both])
        p_solo = len(p_xor_dp[p_xor_dp])
        n_solo = len(n_xor_dn[n_xor_dn])
        p_none = len(p_neither[p_neither])
        n_none = len(n_neither[n_neither])
        print('N sims: ' + str(len(self.results)))
        print('Proton both:    ' + str(p_both) + ', ' + str(p_both / len(self.results) * 100.) + '%')
        print('One, not both:  ' + str(p_solo) + ', ' + str(p_solo / len(self.results) * 100.) + '%')
        print('Proton none:    ' + str(p_none) + ', ' + str(p_none / len(self.results) * 100.) + '%')
        print()
        print('Neutron both:   ' + str(n_both) + ', ' + str(n_both / len(self.results) * 100.) + '%')
        print('One, not both:  ' + str(n_solo) + ', ' + str(n_solo / len(self.results) * 100.) + '%')
        print('Neutron none:   ' + str(n_none) + ', ' + str(n_none / len(self.results) * 100.) + '%')

        print()
        print(np.asarray(p_list)[p_xor_dp])
        print()
        print(np.asarray(dp_list)[p_xor_dp])
        print()
        print(np.asarray(n_list)[n_xor_dn])
        print()
        print(np.asarray(dn_list)[n_xor_dn])
        print()

        if (len(np.asarray(self.results)[p_neither])>0 or len(np.asarray(self.results)[n_neither])>0):                
            fig0 = plt.figure(figsize=[15,15])
            print()
            print('Neither (proton): ')
            for r in np.asarray(self.results)[p_neither]:
                p = r.getEarthRadii(r.p_telemetry[-1])
                pl = r.getEarthRadii(r.p_last)
                dp = r.getEarthRadii(r.dp_telemetry[-1])
                dpl = r.getEarthRadii(r.dp_last)
                print('\t' + r.filename + ': ')
                print('\t\t' + str(p)  + ' => ' + str(pl))
                print('\t\t' + str(dp) + ' => ' + str(dpl))
                pos_p = []
                pos_dp = []
                for t in zip(r.p_telemetry, r.dp_telemetry):
                    pos_p.append(r.getEarthRadii(t[0]))
                    pos_dp.append(r.getEarthRadii(t[1]))
                plt.plot(pos_p)
                plt.plot(pos_dp)
            print()
            print('Neither (neutron): ')
            for r in np.asarray(self.results)[n_neither]:
                n = r.getEarthRadii(r.n_telemetry[-1])
                nl = r.getEarthRadii(r.n_last)
                dn = r.getEarthRadii(r.dn_telemetry[-1])
                dnl = r.getEarthRadii(r.dn_last)
                print('\t' + r.filename + ': ')
                print('\t\t' + str(n)  + ' => ' + str(nl))
                print('\t\t' + str(dn) + ' => ' + str(dnl))
                pos_n = []
                pos_dn = []
                for t in zip(r.n_telemetry, r.dn_telemetry):
                    pos_n.append(r.getEarthRadii(t[0]))
                    pos_dn.append(r.getEarthRadii(t[1]))
                plt.plot(pos_n)
                plt.plot(pos_dn)
            plt.xlim(.9,1.1)
        
        #fig1 = plt.figure(figsize=[15,15])
        #bins = np.linspace(0., 1 + 100 * atol, 100)
        #plt.hist(p_list,  bins=bins, density=True, log=True, color=mpl.colors.to_rgba('b',.3), label='Proton')
        #plt.hist(dp_list, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('m',.3), label='Z-1 Daughter')
        #plt.hist(n_list,  bins=bins, density=True, log=True, color=mpl.colors.to_rgba('r',.3), label='Neutron')
        #plt.hist(dn_list, bins=bins, density=True, log=True, color=mpl.colors.to_rgba('y',.3), label='Z Daughter')
        #plt.legend()
        #plt.show()
        
        p_dp_dist = []
        for pair in np.asarray(self.results)[p_dp_both]:
            p_dp_dist.append(Results.HaversineSeparation(pair.p_last, pair.dp_last))
            
        n_dn_dist = []
        for pair in np.asarray(self.results)[n_dn_both]:
            n_dn_dist.append(Results.HaversineSeparation(pair.n_last, pair.dn_last))
        
        p_dp_dist = np.asarray(p_dp_dist) * units.SI.radius_earth
        n_dn_dist = np.asarray(n_dn_dist) * units.SI.radius_earth
        
        fig2 = plt.figure(figsize=[15,15])
        dist = np.concatenate([p_dp_dist, n_dn_dist])
        n_under1 = len(dist[dist<1.])
        n_under10 = len(dist[dist<10.])
        n_under50 = len(dist[dist<50.])        
        print()
        print('Fraction under 1 m:  ' + str(n_under1/len(dist)*100.))
        print('Fraction under 10 m: ' + str(n_under10/len(dist)*100.))
        print('Fraction under 50 m: ' + str(n_under50/len(dist)*100.))
        non_zero = dist[dist > 0.]
        lo_x = np.log10(min(non_zero) / 100.)
        np.place(dist, dist==0., lo_x)
        hi_x = np.log10(max(dist) * 10.)
        bins = np.logspace(lo_x, hi_x, 100)
        n1, b, p = plt.hist(p_dp_dist, bins=bins, density=False, log=True, color=mpl.colors.to_rgba('b',.3), label='proton-daughter')
        n2, b, p = plt.hist(n_dn_dist, bins=bins, density=False, log=True, color=mpl.colors.to_rgba('r',.3), label='neutron-daughter')
        n = np.concatenate((n1[n1>0], n2[n2>0]))
        lo = np.min(n) / 2.
        hi = np.max(n) * 2.
        plt.plot([1e0,1e0],[lo, hi],'k')
        plt.plot([1e1,1e1],[lo, hi],'k--')
        plt.plot([1e2,1e2],[lo, hi],'k.')
        plt.xscale('log')
        plt.xlim(bins[0], bins[-1])
        plt.ylim(lo, hi)
        plt.legend()
        plt.show()