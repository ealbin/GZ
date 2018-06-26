#!/usr/bin/env python

"""
1) Create masterlist of jobs: Jobs.masterlist()
    (sweep parameters are located here)
2) Write individual job files: Jobs.bundle()
3) Kick-off jobs on the gpatlas cluster: Jobs.submit()
"""

import numpy as np
import os
import subprocess
import sys
import time

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"

################################################
# Edit this:
data_directory = '/beegfs/DATA/atlas/ealbin/GZ'
################################################

jobs_directory = os.path.join(data_directory, 'jobs')
if not os.path.isdir(jobs_directory):
    os.makedirs(jobs_directory)
masterlistpath = os.path.join(jobs_directory, 'master_list.txt')

def masterlist(N_nodes=7, OVERWRITE=False):
    """Create masterlist of jobs to be done.
    N_nodes is only used for estimating time to completion.
    OVERWRITE = True will overwrite existing computations.
    """
    global data_directory, masterlistpath
    
    N_nodes = 7
    N_cpus_per_node = 16
    N_cpus = N_nodes * N_cpus_per_node
    OVERWRITE = False
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

    ##################################
    #       Sweep parameters
    #---------------------------------
    Z_list = np.array([Z_neutron, Z_proton, Z_alpha, Z_oxygen, Z_iron, Z_uranium])

    energy_list = np.array([2e18, 20e18, 200e18])

    radii_list  = np.linspace(0, 6, 61)[::5][1:-1]      # 0.5, 1.0, ..., 5.0, 5.5 [AU]
    theta_list  = np.linspace(0, np.pi, 21)[::4][1:-1]  # 36, 72, 108, 144 [deg] (actually radians)
    phi_list    = np.linspace(0, 2*np.pi, 21)[::5][:-1] # 0, 90, 180, 270 [deg] (actually radians)
    #---------------------------------
    ##################################

    #N_combinations = ( Z_list.size * energy_list.size * 
    #                   radii_list.size * theta_list.size * phi_list.size )
    #print 'N combinations = {0}'.format(N_combinations)

    #est_filesize = 0.5      # [MB]
    #est_elapsed_time = 240. # [seconds]
    #print 'Estimated completed filesize: {0} [GB]'.format(est_filesize * N_combinations / 1000.)
    #print 'Estimated time to completion: {0} [hrs]'.format(N_combinations / float(N_cpus) * est_elapsed_time / 3600.)

    # organize by:
    # GZ/energy/radius/theta_phi_Z.dat
    with open(masterlistpath,'w') as f:
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
                            if os.path.exists(filepath) and (not OVERWRITE):
                                continue
                            start_x = R * np.sin(T) * np.cos(P)
                            start_y = R * np.sin(T) * np.sin(P)
                            start_z = R * np.cos(T)
                    
                            f.write( '{0} {1} {2} {3} {4} {5}\n'.format(start_x, start_y, start_z, Z, E, filepath) )
                        
def bundle(N_jobs=1000, N_paths_per_job=10): 
    """Write individual job files for submitting to gpatlas.
    N_jobs is the maximum number of job files that will be prepared.
    N_paths_per_job is the number of trajectories calculated per bulk job.
    """
    global jobs_directory, masterlistpath

    with open(masterlistpath) as master_file:
        masterlist = master_file.readlines()
    masterlist.reverse()
    
    jobs_submitted = 0    
    while (jobs_submitted < N_jobs) and (len(masterlist) > 0):
        jobs_submitted += 1
        
        with open( os.path.join(jobs_directory, 'job_{0:06}.txt'.format(jobs_submitted)), 'w' ) as job_file:
            for i in xrange(N_paths_per_job):
                if len(masterlist) == 0:
                    break
                job_file.write( masterlist.pop() )

    masterlist.reverse()
    with open(masterlistpath, 'w') as master_file:
        for job in masterlist:
            master_file.write(job)
    

def submit(job_cap=None):
    """Submit jobs to slurm for gpatlas processing.
    If job_cap=None, all jobs existing will be submitted.
    Setting job_cap = [integer] will submit the first [integer] number of jobs.
    """
    global jobs_directory
    
    jobs = os.listdir(jobs_directory)
    jobs.sort()
    jobs = [j for j in jobs if j.find('job_') == 0]
    
    for i, job in enumerate(jobs):
        if job_cap is not None:
            if i+1 > job_cap:
                break
        subprocess.call(['sbatch', '../bin/slurm_batch.sh', os.path.join(jobs_directory, job) ])
        time.sleep(0.1)
