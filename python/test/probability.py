#!/user/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

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