#!/user/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time

import gz

def run(position=[1,0,0], beta=[-1,0,0], E=2e18, samples=100000):
    
    Zlist = [2, 26, 92]
    
    colors=[(1,0,0), (0,1,0), (0,0,1)]
    
    plt.figure(figsize=[15,15])
    plt.xscale('log')
    ebins = np.logspace(-3, 2, 100)
    for _, z in enumerate(Zlist):
        gammas, ee, ii = gz.probability.Solar.get_photon(position, beta, z, E, size=samples, plottables=True)
        label = 'Z=' + str(z) + ', ' + 'E=' + str(E)
        plt.hist(gammas, bins=ebins, density=True, range=[1e-4, 10], log=True, color=colors[_] + (.3,), label=label)
        plt.plot(ee, ii, color=colors[_], linewidth=2)
    plt.legend(fontsize='xx-large')
    plt.xlim(ebins[0], ebins[-1])
    plt.ylim(1e-6, 1e1)
    
    
def longshot(position=[1,0,0], beta=[1,0,0], Z=92, E=1e18, R_limit=40, outname='out_test'):

    outpath = gz.earth.Earth.OUT_JOB_PATH + gz.path.Outgoing.DEFAULT_SAVE_PATH[1:]
    subdir = str(np.abs(Z)) + '_' + str(int(E/1e15))
    
    print('Loading HMF interpolation maps... ', end='', flush=True)
    start = time.time()
    _ = gz.magnetic_field.cartesianTesla(gz.coordinates.Cartesian.earth)
    print('{:.5f} sec'.format(time.time() - start))
    
    print('Working on outgoing... ', end='', flush=True)
    start = time.time()
    outgoing = gz.path.Outgoing(position, beta, Z, E, max_step=.01, 
                                R_limit=R_limit, save_path=outpath, filename=outname)
    outgoing.propagate()
    print('{:.5f} sec'.format(time.time() - start))
    
    outfile =  outpath + '/' + subdir + '/' + outname + '.outgoing'
    
    print('Working on probability... ', flush=True)
    start = time.time()
    gz.earth.Earth.incoming_jobs(filelist=[outfile], runs=1, plot=True, histograms=False)
    print('time: {:.5f} sec'.format(time.time() - start))
    
    print('finished')
    print("(don't forget to delete output file before re-running)")
    print(outfile)