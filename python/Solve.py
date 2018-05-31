import numpy as np
import platform
from scipy import integrate
import sys
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import Dynamics

# physical constants / conversions
m_per_AU = 149597870700. # unit conversion [meters / astronomical unit]
Rs   = 695508000.        # radius of the Sun [meters]
RsAU = Rs / m_per_AU     # radius of the Sun [AU]
Re   = 6378100.          # radius of Earth [meters]
ReAU = Re / m_per_AU     # radius of Earth [AU]
c    = 299792458.        # speed of light [meters / second]

def trajectory( start_pos, Z, E, savefile='./path.data', start_beta=None, 
                algorithm='vode', stepsize=100000):
    """Calculates the trajectory of nuclei Z with energy E as it travels 
    throught the solar system. Assumes ultra-relativistic nuclei, aka
    |beta| > 0.999990 ~= 1. 
        start_pos === start position in solar cartesian coordinates [astronomical units]
        Z === number of protons in nuclei [unitless]
        E === energy of nuclei [electronVolts] === (2..200)e18 eV
        savefile === textfile to save trajectory positions
        start_beta === relativistic velocity [unit-less] (default is towards Earth)
        algorithm === default "vode" (best performance), "lsoda" (almsot equal performance),
                      "dop853" (worst performance, explicit runge-kutta).
        stepsize === integration stepsize [meters] (default 100 km)
    """
    start_time = time.time()
    initial_pos = np.array(start_pos)

    # figure out beta
    initial_beta = None
    if start_beta is None:
        cartesian_earth = np.asarray([1,0,0])
        initial_beta = cartesian_earth - initial_pos # (un-normalized heading)
    else:
        initial_beta = np.asarray(start_beta)
    # enforce special relativity (|beta| ~= 1):    
    initial_beta = initial_beta / np.sqrt( np.dot(initial_beta, initial_beta) )

    initial_conditions = np.concatenate(( initial_pos, initial_beta ))
    dynamic_ratio = Z / float(E)
    
    # configure integrator
    integrator = None
    success_code = None
    if algorithm == 'vode':
        integrator = integrate.ode(Dynamics.applyForces).set_integrator('vode', method='BDF')    
        success_code = 2
    elif algorithm == 'lsoda':
        integrator = integrate.ode(Dynamics.applyForces).set_integrator('lsoda', method='BDF')    
        success_code = 2
    elif algorithm == 'dop853':
        integrator = integrate.ode(Dynamics.applyForces).set_integrator('dop853')    
        success_code = 1
    else:
        print "invalid algorithm, please input 'vode', 'lsoda' or 'dop853'"
        sys.exit(0)

    initial_l = 0. # path begins at 0 distance [AU]
    integrator.set_initial_value(initial_conditions, initial_l).set_f_params(dynamic_ratio)
    dl = stepsize / m_per_AU  # numerical stepsize, dl distance (in AU)
    positions = [] # container for positions computed along the path
    
    def isOutOfBounds(position, spacelimit=6, sun_pos=np.array([0,0,0])):
        """ Checks if computed path has wondered out-of-bounds.
        returns True beyond spacelimit, False if not.
        """
        sun_diff = position - sun_pos
        Rsun     = np.sqrt( np.dot(sun_diff, sun_diff) )
        if Rsun > spacelimit:
            return True
        return False
        
    def isNearSun(position, sunlimit=2*RsAU, sun_pos=np.array([0,0,0])):
        """ Checks if computed path has found its way to the Sun.
        returns True if near the Sun, False if not.
        """
        sun_diff = position - sun_pos
        Rsun     = np.sqrt( np.dot(sun_diff, sun_diff) )
        if Rsun < sunlimit:
            return True
        return False
        
    def isNearEarth(position, earthlimit=2*ReAU, earth_pos=np.array([1,0,0])):
        """ Checks if computed path has found its way to Earth.
        returns True if near Earth, False if not.
        """
        earth_diff = position - earth_pos
        Rearth     = np.sqrt( np.dot(earth_diff, earth_diff) )
        if Rearth < earthlimit:
            return True
        return False
    
    exit_status = None
    while integrator.successful():
        integrator.integrate( integrator.t + dl ) # because "dt" is "dl" in this context
        positions.append( np.asarray(integrator.y[:3]) ) # computed position [AU]
        if isOutOfBounds(positions[-1]):
            exit_status = 'out-of-bounds'
            break
        if isNearSun(positions[-1]):
            exit_status = 'near-sun'
            break
        if isNearEarth(positions[-1]):
            exit_status = 'near-earth'
            break
    
    # thin out data
    sizecap  = 500
    skipsize = 1000
    if len(positions) > 3*sizecap:
        front  = positions[:sizecap]
        middle = positions[sizecap+1:-sizecap-1:skipsize]
        back   = positions[-sizecap:]
        positions = front + middle + back
    
    positions = np.asarray(positions)
    x, y, z = positions.T

    # save to disk
    with open(savefile, 'w') as file:
        file.write('whatami: result from Solve.trajectory\n')
        file.write('system_info: {}\n'.format(platform.uname()))
        file.write('elapsed_time[sec]: {}\n'.format(time.time()-start_time))
        file.write('algorithm: {}\n'.format(algorithm))
        file.write('stepsize[m]: {}\n'.format(stepsize))
        if integrator.get_return_code() == success_code:
            file.write('integration_status: successful\n')
        else:
            file.write('integration_status: ***failed***\n')
        file.write('exit_status: {}\n'.format(exit_status))
        file.write('number_of_protons_Z: {}\n'.format(Z))
        file.write('energy[eV]: {}\n'.format(E))
        file.write('start_pos[AU]: ')
        np.savetxt(file, initial_pos, delimiter=',', newline=',')
        file.write('\n')
        file.write('start_beta[unit-less]: ')
        np.savetxt(file, initial_beta, delimiter=',', newline=',')
        file.write('\n')
        file.write('x_positions[AU]: ')
        np.savetxt(file, x, delimiter=',', newline=',')
        file.write('\n')
        file.write('y_positions[AU]: ')
        np.savetxt(file, y, delimiter=',', newline=',')
        file.write('\n')
        file.write('z_positions[AU]: ')
        np.savetxt(file, z, delimiter=',', newline=',')
        file.write('\n')
        file.write('END OF LINE\n')

def plot(filelist=None):
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