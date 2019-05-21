#!/user/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import time

import gz

def run(Z=92, E=2e18, beta=[0, 1, 0], plot=True, positions=None, betas=None, zigzag=False):
    start = time.time()
    
    if (positions is None):
        e_position = gz.coordinates.Cartesian.earth
        d_position = gz.coordinates.Cartesian.earth
    else:
        e_position = positions[:3]
        d_position = positions[3:]
    
    if (betas is None):
        e_beta = beta
        d_beta = beta
    else:
        e_beta = betas[:3]
        d_beta = betas[3:]

    eul_out = gz.path.Outgoing(e_position, e_beta, Z, E, zigzag=zigzag, save=False, max_step=.1)    
    dop_out = gz.path.Outgoing(d_position, d_beta, Z, E, zigzag=zigzag, save=False, max_step=.1)
    
    eul_start = time.time()
    eul_out.propagate(B_override=[0,0,1e-4], algorithm='euler')
    eul_elapsed = time.time() - eul_start
    
    dop_start = time.time()
    dop_out.propagate(B_override=[0,0,1e-4], algorithm='dop853')
    dop_elapsed = time.time() - dop_start

    ex = [t[0] for t in eul_out.telemetry]
    ey = [t[1] for t in eul_out.telemetry]
    
    dx = [t[0] for t in dop_out.telemetry]
    dy = [t[1] for t in dop_out.telemetry]

    if (plot):    
        plt.figure(figsize=[8,8])
        plt.plot(ex, ey, 'b')
        plt.plot(dx, dy, 'g')
        plt.show()
        
        print('Euler elapsed:  ' + str(eul_elapsed))
        print('Dop853 elapsed: ' + str(dop_elapsed))
        print()
        
        difference = eul_out.position - dop_out.position
        difference = np.sqrt( np.dot(difference, difference) )
        print('Difference [AU]:     ' + str(difference))
        difference *= gz.units.Change.AU_to_meter
        print('Difference [meters]: ' + str(difference))
        print()
        
        print('Euler steps:  ' + str(len(eul_out.telemetry)))
        print('Dop853 steps: ' + str(len(dop_out.telemetry)))
    else:
        return ex, ey, eul_out.position, eul_out.beta, dx, dy, dop_out.position, dop_out.beta


def zigzag(Z=92, E=2e18, beta=[0,1,0], trips=10):
    """dosen't work without editing path.py, line 174 to:
        while (self.distance < self.R_limit + Outgoing.LIMIT_BUFFER): # and self.dist_sun < self.R_limit):            
    """
    beta = np.asarray(beta, dtype=np.float64)

    positions = None
    betas = np.concatenate([beta, beta])
    
    ex = []
    ey = []
    dx = []
    dy = []
    for _ in range(2 * trips):
        print('running ' + str(_+1) + ' of ' + str(2*trips))
        exx, eyy, epos, eb, dxx, dyy, dpos, db = run(Z=Z, E=E, betas=betas, plot=False, positions=positions, zigzag=True)
        ex.append(exx)
        ey.append(eyy)
        dx.append(dxx)
        dy.append(dyy)
        Z = -1 * Z
        e_beta = -1. * eb
        d_beta = -1. * db
        positions = np.concatenate([epos, dpos])
        betas = np.concatenate([e_beta, d_beta])

    ex = np.concatenate(ex)
    ey = np.concatenate(ey)
    dx = np.concatenate(dx)
    dy = np.concatenate(dy)
    
    plt.figure(figsize=[8,8])
    plt.plot(ex, ey, 'b')
    plt.plot(dx, dy, 'g')
    plt.show()
    
    earth = gz.coordinates.Cartesian.earth
    e_err = epos - earth
    d_err = dpos - earth
    
    e_err = np.sqrt(np.dot(e_err, e_err))
    d_err = np.sqrt(np.dot(d_err, d_err))
    
    print('Euler error [AU]:  ' + str(e_err))
    print('Euler error [m]:   ' + str(gz.units.Change.AU_to_meter*e_err))
    print()
    print('Dop853 error [AU]: ' + str(d_err))
    print('Dop853 error [m]:  ' + str(gz.units.Change.AU_to_meter*d_err))
    
    
def plot(filelist):
    fig = plt.figure(figsize=[8,8])
    ax = plt.axes(projection='3d')
    phi = None
    
    for file in filelist:
        theta = float(os.path.split(file)[1].split('_')[0])
        phi = float(os.path.split(file)[1].split('_')[1])
        with open(file) as f:
            x = []
            y = []
            z = []

            for line in f.readlines():
                if (line.startswith('#')):
                    continue
                pos_x, pos_y, pos_z, beta_x, beta_y, beta_z, dist = line.split()
                x.append(float(pos_x))
                y.append(float(pos_y))
                z.append(float(pos_z))
            ax.plot3D(x, y, z, 'gray')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.view_init(azim=phi-90, elev=0)
    
