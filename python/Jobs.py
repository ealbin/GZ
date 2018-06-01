"""Kick-off jobs on the gpatlas cluster
"""

import numpy as np
import os
import subprocess
import sys

# cluster information (gpatlas)
data_directory = '/DFS-L/DATA/atlas/ealbin/GZ'
N_nodes = 7
N_cpus_per_node = 16
N_cpus = N_nodes * N_cpus_per_node
print 'N cluster cpus: {0}'.format(N_cpus)

# convienient constants
Z_neutron  = 0
Z_proton   = 1
Z_deuteron = 1
Z_alpha    = 2
Z_helium   = 4
Z_carbon   = 6
Z_nitrogen = 7
Z_oxygen   = 8
Z_neon     = 10
Z_argon    = 18
Z_iron     = 26
Z_krypton  = 36
Z_xenon    = 54
Z_uranium  = 92

# Sweep parameters
#--------------------------------
#--------------------------------
Z_list = np.array([Z_neutron, Z_proton, Z_alpha, Z_oxygen, Z_iron, Z_uranium])

energy_list = np.array([2e18, 20e18, 200e18])

radii_list  = np.linspace(0, 6, 61)[1:-1]      # 0.1, 0.2, ..., 5.8, 5.9 [AU]
theta_list  = np.linspace(0, np.pi, 21)[1:-1]  # 9, 18, 27, ..., 162, 171 [deg] (actually radians)
phi_list    = np.linspace(0, 2*np.pi, 21)[:-1] # 0, 18, 36, ..., 324, 342 [deg] (actually radians)
#--------------------------------
#--------------------------------

N_combinations = ( species_list.size * energy_list.size * 
                   radii_list.size * theta_list.size * phi_list.size )
print 'N combinations = {0}'.format(N_combinations)

est_filesize = 1.       # [MB]
est_elapsed_time = 180. # [seconds]
print 'Estimated completed filesize: {0} [GB]'.format(est_filesize * N_combinations / 1000.)
print 'Estimated time to completion: {0} [hrs]'.format(N_combinations / float(N_cpus) * est_elapsed_time / 3600.)

# organize by:
# GZ/energy/radius/theta_phi_Z.dat
for E in energy_list:
    for R in radii_list:
        dirname  = os.path.join( 'energy_{0}e18'.format(int(E/1e18)), 'radius_{0}'.format(R) )
        fullpath = os.path.join(data_directory, dirname)
        if not os.path.isdir(fullpath):
            os.makedirs(fullpath)

        for T in theta_list:
            for P in phi_list:
                for Z in Z_list:
                    filename = 'theta_{0}_phi_{1}_Z_{2}.dat'.format( int(T*180./np.pi), int(P*180./np.pi), Z )
                    filepath = os.path.join( fullpath, filename )
                    start_x = R * np.sin(T) * np.cos(P)
                    start_y = R * np.sin(T) * np.sin(P)
                    start_z = R * np.cos(T)
                    
                    subprocess.call(['sbatch', 'sbatch_interface.sh', start_x, start_y, start_z, Z, E, filepath])
                    sys.exit()