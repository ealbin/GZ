#!/usr/bin/env python

"""A collection of useful scripts to post-process simulation results.
"""

import numpy as np


__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"


def read(file):
    """Reads a data file and returns a dictionary of its data.
    """
    file_dict = {'name':file}
    with open(file) as f:
        # chew first line (header)
        line = f.readline()

        # system info
        line = f.readline()
        system_info = line.split(':')[-1].strip()
        file_dict['system'] = system_info

        # elapsed time [seconds]
        line = f.readline()
        elapsed = float( line.split(':')[-1].strip() )
        file_dict['time'] = elapsed

        # integration algorithm
        line = f.readline()
        algorithm = line.split(':')[-1].strip()
        file_dict['algorithm'] = algorithm

        # stepsize [meters]
        line = f.readline()
        stepsize = float( line.split(':')[-1].strip() )
        file_dict['stepsize'] = stepsize

        # integration status
        line = f.readline()
        int_status = line.split(':')[-1].strip()
        file_dict['integration'] = int_status

        # exit status
        line = f.readline()
        exit_status = line.split(':')[-1].strip()
        file_dict['exit_info'] = exit_status

        # distance from earth
        line = f.readline()
        earth_dist = float( line.split(':')[-1].strip() )
        file_dict['earth_dist'] = earth_dist
        
        # number of protons (Z)
        line = f.readline()
        Z = float( line.split(':')[-1].strip() )
        file_dict['Z'] = Z

        # energy (E) [electronVolts]
        line = f.readline()
        E = float( line.split(':')[-1].strip() )
        file_dict['E'] = E
        
        # start position [AU]
        line = f.readline()
        start_pos = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])
        file_dict['initial_pos'] = start_pos
        
        # start beta [unit-less]
        line = f.readline()
        start_beta = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])
        file_dict['initial_beta'] = start_beta
        
        # x positions [AU]
        line = f.readline()
        pos_x = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])
        file_dict['x'] = pos_x
        
        # y positions [AU]
        line = f.readline()
        pos_y = np.asarray([ float(y) for y in line.split(':')[-1].strip().strip(',').split(',') ])
        file_dict['y'] = pos_y
        
        # z positions [AU]
        line = f.readline()
        pos_z = np.asarray([ float(z) for z in line.split(':')[-1].strip().strip(',').split(',') ])
        file_dict['z'] = pos_z
        
        # final position [AU]
        stop_pos = np.asarray([ pos_x[-1], pos_y[-1], pos_z[-1] ])
        file_dict['last_pos'] = stop_pos
        
        # final beta [unit-less]
        stop_beta = stop_pos - np.asarray([ pos_x[-2], pos_y[-2], pos_z[-2] ])
        stop_beta = stop_beta / np.sqrt( np.dot(stop_beta, stop_beta) )    
        file_dict['last_beta'] = stop_beta
    
        # heading change [degrees]
        delta_heading = np.arccos(np.dot( start_beta, stop_beta )) * 180. / np.pi
        file_dict['delta_heading'] = delta_heading

    return file_dict
