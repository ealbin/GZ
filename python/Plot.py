
import numpy as np
import os

def select(data_path, energy_range=None, radius_range=None, theta_range=None, phi_range=None, Z_range=None):
    """Collect trajectories matching criteria.
    data_path: base directory for trajectory data.
    *_range: if None, select all.  Otherwise select between either [low, high], 
             or [[val1, val2, ..., valN]] specifically.
    e.g. select radius=0.1, 1.1, 2.1, theta=45 to 135, and all energy, phi and Z:
         radius_range=[[0.1, 1.1, 2.1]], theta_range=[45, 135]
    returns a list of files satisfying the requested criteria.
    """
    def extractNumber(substring):
        slash_pos = substring.find('/')
        under_pos = substring.find('_')
        stop_pos  = min([slash_pos, under_pos])
        return float(substring[:stop_pos])
    
    def inRange(val, val_range):
        val_range = np.asarray(val_range)
        shape     = val_range.shape
        if shape[0] == 1:
            # pick exact values
            if val in val_range:
                return True
        else:
            # specified range
            if (val_range.min() < val) and (val < val_range.max()):
                return True
        return False
        
    energy_str = 'energy_'
    radius_str = 'radius_'
    theta_str  = 'theta_'
    phi_str    = 'phi_'
    Z_str      = 'Z_'
    
    file_list = []
    for root, dirs, files in os.walk(data_path):
        for name in files:
            file       = os.path.join(root, name)
            energy_pos = file.find(energy_str) + len(energy_str)
            radius_pos = file.find(radius_str) + len(radius_str)
            theta_pos  = file.find(theta_str)  + len(theta_str)
            phi_pos    = file.find(phi_str)    + len(phi_str)
            Z_pos      = file.find(Z_str)      + len(Z_str)
    
            # skip non-data files
            if ( (energy_pos == -1) or (radius_pos == -1) or
                 (theta_pos  == -1) or (phi_pos    == -1) or (Z_pos == -1) ):
                continue
            
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
            
    return file_list

def plot(file_list, xlim=None, ylim=None, zlim=None):
    # these are not available on gpatlas
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig_full = plt.figure(figsize=(15,10))
    ax_full  = fig_full.gca(projection='3d')
    ax_full.set_xlim(-6,6)
    ax_full.set_ylim(-6,6)
    ax_full.set_zlim(-6,6)
    ax_full.set_xlabel('x [AU]')
    ax_full.set_ylabel('y [AU]')
    ax_full.set_zlabel('z [AU]')
    ax_full.set_title('Trajectories')

    fig_earth = plt.figure(figsize=(15,10))
    ax_earth  = fig_earth.gca(projection='3d')
    ax_earth.set_xlim(-10,10) # factors of Earth Radius
    ax_earth.set_ylim(-10,10)
    ax_earth.set_zlim(-10,10)
    ax_earth.set_xlabel('x [Earth Radii]')
    ax_earth.set_ylabel('y [Earth Radii]')
    ax_earth.set_zlabel('z [Earth Radii]')
    ax_earth.set_title('Near Earth Trajectories')

    # draw earth
    phi, theta = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
    x = ReAU * np.sin(theta) * np.cos(phi) + 1.
    y = ReAU * np.sin(theta) * np.sin(phi)
    z = ReAU * np.cos(theta)
    ax_full.plot_wireframe(x, y, z, color="b")    
    ax_earth.plot_wireframe((x-1.)/ReAU, y/ReAU, z/ReAU, color="b")    
    
    fig_sun = plt.figure(figsize=(15,10))
    ax_sun  = fig_sun.gca(projection='3d')
    ax_sun.set_xlim(-50, 50)  # factors of Sun Radius
    ax_sun.set_ylim(-50, 50)
    ax_sun.set_zlim(-50, 50)
    ax_sun.set_xlabel('x [Sun Radii]')
    ax_sun.set_ylabel('y [Sun Radii]')
    ax_sun.set_zlabel('z [Sun Radii]')
    ax_sun.set_title('Near Sun Trajectories')

    # draw sun
    phi, theta = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = RsAU * np.sin(theta) * np.cos(phi)
    y = RsAU * np.sin(theta) * np.sin(phi)
    z = RsAU * np.cos(theta)
    ax_full.plot_wireframe(x, y, z, color="y")    
    ax_sun.plot_wireframe(x/RsAU, y/RsAU, z/RsAU, color="y")    
    
    for file in filelist:
        with open(file) as f:
    
            # chew first line (header)
            line = f.readline()

            # system info
            line = f.readline()
            system_info = line.split(':')[-1].strip()

            # elapsed time [seconds]
            line = f.readline()
            elapsed = float( line.split(':')[-1].strip() )

            # integration algorithm
            line = f.readline()
            algorithm = line.split(':')[-1].strip()

            # stepsize [meters]
            line = f.readline()
            stepsize = float( line.split(':')[-1].strip() )

            # integration status
            line = f.readline()
            int_status = line.split(':')[-1].strip()
            # check if == 'successful'

            # exit status
            line = f.readline()
            exit_status = line.split(':')[-1].strip()
            # check if == 'near-earth'

            # number of protons (Z)
            line = f.readline()
            Z = int( line.split(':')[-1].strip() )

            # energy (E) [electronVolts]
            line = f.readline()
            E = float( line.split(':')[-1].strip() )

            # start position [AU]
            line = f.readline()
            start_pos = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])

            # start beta [unit-less]
            line = f.readline()
            start_beta = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])

            # x positions [AU]
            line = f.readline()
            pos_x = np.asarray([ float(x) for x in line.split(':')[-1].strip().strip(',').split(',') ])

            # y positions [AU]
            line = f.readline()
            pos_y = np.asarray([ float(y) for y in line.split(':')[-1].strip().strip(',').split(',') ])

            # z positions [AU]
            line = f.readline()
            pos_z = np.asarray([ float(z) for z in line.split(':')[-1].strip().strip(',').split(',') ])

            # final position [AU]
            stop_pos = np.asarray([ pos_x[-1], pos_y[-1], pos_z[-1] ])

            # final beta [unit-less]
            stop_beta = stop_pos - np.asarray([ pos_x[-2], pos_y[-2], pos_z[-2] ])
            stop_beta = stop_beta / np.sqrt( np.dot(stop_beta, stop_beta) )    
    
            # heading change [degrees]
            delta_heading = np.arccos(np.dot( start_beta, stop_beta )) * 180. / np.pi

            print 'file: {}'.format(file)
            print 'elapsed time: {}'.format(elapsed)
            print 'exit statii: {}, {}'.format(int_status, exit_status)
            print 'Z, E: {}, {}'.format(Z, E)
            last_beta = np.array([pos_x[-1], pos_y[-1], pos_z[-1]]) - np.array([pos_x[-2], pos_y[-2], pos_z[-2]])
            last_beta = last_beta / np.sqrt(np.dot(last_beta, last_beta))
            print 'Heading change [degrees]: {}'.format( np.arccos(start_beta, last_beta) )
            last_pos_earth = np.array([pos_x[-1] - 1., pos_y[-1], pos_z[-1]]) 
            print 'End distance from Earth [Earth Radii]: {}'.format( np.sqrt(np.dot(last_pos_earth, last_pos_earth))/ReAU )
            print 'End coordinate location [AU]: {}, {}, {}'.format( pos_x[-1], pos_y[-1], pos_z[-1] )
            print 'End coordinate location from Earth [Earth Radii]: {}'.format(last_pos_earth/ReAU)            
            print
            
            ax_full.plot(pos_x, pos_y, pos_z)
            ax_earth.plot((pos_x-1.)/ReAU, pos_y/ReAU, pos_z/ReAU)
            ax_sun.plot(pos_x/RsAU, pos_y/RsAU, pos_z/RsAU)

    plt.show()        