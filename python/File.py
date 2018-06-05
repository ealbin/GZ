
import numpy as np
import os
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# physical constants / conversions
m_per_AU = 149597870700. # unit conversion [meters / astronomical unit]
Rs   = 695508000.        # radius of the Sun [meters]
RsAU = Rs / m_per_AU     # radius of the Sun [AU]
Re   = 6378100.          # radius of Earth [meters]
ReAU = Re / m_per_AU     # radius of Earth [AU]

def select(data_path, energy_range=None, radius_range=None, 
           theta_range=None, phi_range=None, Z_range=None, special=None):
    """Collect trajectories matching criteria.
    data_path: base directory for trajectory data.
    *_range: if None, select all.  Otherwise select scalar number, or [low, high] (inclusive), 
             or [[val1, val2], [val3,val4]] for multiple ranges.
    e.g. select radius=0.1, theta=45 to 135, phi=0 to 90 and 270 to 360, and all energy and Z:
         radius_range=0.1, theta_range=[45, 135], phi_range=[ [0, 90], [270, 360] ]
    special: automatically filters (over-rules theta/phi ranges if specified) based on
             'front-cone': theta_range=[60, 120]; phi_range=[ [0, 30], [330, 360] ]
             'back-cone' : theta_range=[60, 120]; phi_range=[150, 210]
             'front-lobe': theta_range=[ [0, 60], [120, 180] ]; phi_range=[ [30, 90], [270, 330] ] 
             'back-lobe' : theta_range=[ [0, 60], [120, 180] ]; phi_range=[ [90, 150], [210, 270] ]
    returns a list of files satisfying the requested criteria.
    """
    def extractNumber(substring):
        slash_pos = substring.find('/')
        under_pos = substring.find('_')
        dot_pos   = substring.find('.')
        if slash_pos == -1:
            slash_pos = len(substring)
        if under_pos == -1:
            under_pos = len(substring)
        if dot_pos == -1:
            dot_pos = len(substring)
        stop_pos  = min([slash_pos, under_pos, dot_pos])
        return float(substring[:stop_pos])
    
    def inRange(val, val_range):
        val_range = np.asarray(val_range)
        shape     = val_range.shape
        shape_len = len(shape)
        if shape_len == 0:
            # it's a number
            if val == val_range:
                return True
        elif shape_len == 1:
            # it's a range [low, high]
            if (val_range.min() <= val) and (val <= val_range.max()):
                return True
        elif shape_len == 2:
            # it's multiple ranges [ [low1, high1], ..., [lowN, highN] ]
            for subrange in val_range:
                if (subrange.min() <= val) and (val <= subrange.max()):
                    return True
        return False

    if special is not None:
        if special == 'front-cone':
            theta_range = [60, 120]
            phi_range   = [ [0, 30], [330, 360] ]
        elif special == 'back-cone':
            theta_range = [60, 120]
            phi_range   = [150, 210]
        elif special == 'front-lobe':
            theta_range = [ [0, 60], [120, 180] ]
            phi_range   = [ [30, 90], [270, 330] ]
        elif special == 'back-lobe':
            theta_range = [ [0, 60], [120, 180] ]
            phi_range   = [ [90, 150], [210, 270] ]
        else:
            print "special = {'front-cone', 'back-cone', 'front-lobe', or 'back-lobe'}"
            sys.exit(1)
    
    energy_str = 'energy_'
    radius_str = 'radius_'
    theta_str  = 'theta_'
    phi_str    = 'phi_'
    Z_str      = 'Z_'
    
    file_list = []
    for root, dirs, files in os.walk(data_path):
        for name in files:
            file       = os.path.join(root, name)
            energy_pos = file.find(energy_str)
            radius_pos = file.find(radius_str)
            theta_pos  = file.find(theta_str)
            phi_pos    = file.find(phi_str)
            Z_pos      = file.find(Z_str)
    
            # skip non-data files
            if ( (energy_pos == -1) or (radius_pos == -1) or
                 (theta_pos  == -1) or (phi_pos    == -1) or (Z_pos == -1) ):
                continue
            
            energy_pos += len(energy_str)
            radius_pos += len(radius_str)
            theta_pos  += len(theta_str)
            phi_pos    += len(phi_str)
            Z_pos      += len(Z_str)
            
            energy = extractNumber( file[energy_pos:] )
            radius = extractNumber( file[radius_pos:] )
            theta  = extractNumber( file[theta_pos:] )
            phi    = extractNumber( file[phi_pos:] )
            Z      = extractNumber( file[Z_pos:] )
    
            if energy_range is not None:
                if not inRange(energy, energy_range):
                    continue
            if radius_range is not None:
                if not inRange(radius, radius_range):
                    continue
            if theta_range is not None:
                if not inRange(theta, theta_range):
                    continue
            if phi_range is not None:
                if not inRange(phi, phi_range):
                    continue
            if Z_range is not None:
                if not inRange(Z, Z_range):
                    continue
            
            file_list.append(file)

    file_list.sort()
    return file_list


def read(file):
    """Reads a data file and returns a dictionary of its data.
    """
    file_dict = {}
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

        # number of protons (Z)
        line = f.readline()
        Z = int( line.split(':')[-1].strip() )
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


def separation(file_list, Z1, Z2=0):
    """Find separation distance [meters] between two species Z1 and Z2 (neutron by default).
    Returns a list of these distances.
    """
    global m_per_AU
    
    Z1_files = [ file for file in file_list if file.find('Z_{0}'.format(Z1)) > -1 ]
    Z2_files = [ file for file in file_list if file.find('Z_{0}'.format(Z2)) > -1 ]
    
    dist_list = []
    for file1 in Z1_files:
        match_pos = file1.find('Z_{0}'.format(Z1))
        match_str = file1[:match_pos]
        for file2 in Z2_files:
            test_str = file2[:match_pos]
            if test_str == match_str:
                z1_dict = read(file1)
                if (z1_dict['exit_info'] != 'near-earth') or (z1_dict['integration'] != 'successful'):
                    break
                z1_pos = z1_dict['last_pos']
                del z1_dict
                
                z2_dict = read(file2)
                if (z2_dict['exit_info'] != 'near-earth') or (z2_dict['integration'] != 'successful'):
                    break
                z2_pos = z2_dict['last_pos']
                del z2_dict 
                
                diff = z2_pos - z1_pos
                dist_list.append( np.sqrt(np.dot(diff, diff)) )
                break
    
    return np.asarray(dist_list) * m_per_AU
    
    
def plot(file_list, ax=None, color=None):
    """
    """
    global ReAU, RsAU
    
    if ax is None:
        fig = plt.figure(figsize=(15,10))
        ax  = fig.gca(projection='3d')

    if color is None:
        color = 'k'
        
    # draw earth
    phi, theta = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
    x = ReAU * np.sin(theta) * np.cos(phi) + 1.
    y = ReAU * np.sin(theta) * np.sin(phi)
    z = ReAU * np.cos(theta)
    ax.plot_wireframe(x, y, z, color="b")    
    
    # draw sun
    phi, theta = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = RsAU * np.sin(theta) * np.cos(phi)
    y = RsAU * np.sin(theta) * np.sin(phi)
    z = RsAU * np.cos(theta)
    ax.plot_wireframe(x, y, z, color="y")    
    
    for file in filelist:
        ax.plot(pos_x, pos_y, pos_z, color=color)

    return ax
