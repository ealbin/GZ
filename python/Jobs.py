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
    thetax_list = np.array([5., 30., 60., 90., 120., 150., 175.]) * np.pi / 180. # [radians]

    Z_list = np.array([Z_helium, Z_oxygen, Z_iron, Z_uranium])
    A_list = np.array([      4.,      16.,    56.,      238.])
    
    energy_list = np.array([2e18, 200e18]) # [eV]
    phix_list   = np.array([0., 120., 240.]) * np.pi / 180. # [radians]
    radii_list  = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]) # [AU]
    
    #---------------------------------
    ##################################

    N_combinations = ( Z_list.size * energy_list.size * 
                       radii_list.size * thetax_list.size * phix_list.size )
    print 'N combinations = {0}'.format(N_combinations)

    # organize by:
    # GZ/thetax/Z_energy_phix_radius.dat  nuclide with A-1 (take half proton, half neutron, aka Z - 0.5 and A - 1)
    # GZ/thetax/Z_energy_phix_radius.pdat proton
    jobcount = 0
    with open(masterlistpath,'w') as f:
        for Tx in thetax_list:
            dirname  = 'thetax_{0}'.format(Tx*180./np.pi)
            fullpath = os.path.join(data_directory, dirname)        
            if not os.path.isdir(fullpath):
                os.makedirs(fullpath)

            for Z, A in zip(Z_list, A_list):
                for E in energy_list:
                    for Px in phix_list:
                        for R in radii_list:
                            skip  = False
                            skipp = False
                            
                            filename = 'Z{0:02}_E{1:02}e18_PHIX_{2:03.0f}_R{3:2.1f}.dat'.format( Z, E/2e18, Px*180./np.pi, R )
                            filepath = os.path.join( fullpath, filename )
                            if os.path.exists(filepath) and (not OVERWRITE):
                                skip = True

                            pfilename = 'Z{0:02}_E{1:02}e18_PHIX_{2:03.0f}_R{3:2.1f}.pdat'.format( Z, E/2e18, Px*180./np.pi, R )
                            pfilepath = os.path.join( fullpath, pfilename )
                            if os.path.exists(pfilepath) and (not OVERWRITE):
                                skip = True

                            start_x = R * np.cos(Tx)
                            start_y = R * np.sin(Tx) * np.sin(Px)
                            start_z = R * np.sin(Tx) * np.cos(Px)
                    
                            # !!!! SIMULATE Z - .5 for DECAY, aka ave between loss of proton vs neutron !!!!
                            if not skip:
                                f.write( '{0} {1} {2} {3} {4} {5}\n'.format(start_x, start_y, start_z, Z -.5, E*(A-1.)/A,  filepath) )
                                jobcount += 1
                            if not skipp:
                                f.write( '{0} {1} {2} {3} {4} {5}\n'.format(start_x, start_y, start_z,     1,        E/A, pfilepath) )
                                jobcount += 1
    print 'prepared {0} jobs'.format(jobcount)
                        
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
