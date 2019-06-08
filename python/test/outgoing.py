#!/user/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import time

import gz

def run(Z=92, E=2e18, beta=[0, 1, 0], B_override=[0,0,1e-4], plot=True, positions=None, betas=None, zigzag=False, interpolate_B=True):
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

    eul_out = gz.path.Outgoing(e_position, e_beta, Z, E, zigzag=zigzag, save=False, max_step=.01)    
    dop_out = gz.path.Outgoing(d_position, d_beta, Z, E, zigzag=zigzag, save=False, max_step=.01)
    
    eul_start = time.time()
    eul_out.propagate(B_override=B_override, interpolate_B=interpolate_B, algorithm='euler')
    eul_elapsed = time.time() - eul_start
    
    dop_start = time.time()
    dop_out.propagate(B_override=B_override, interpolate_B=interpolate_B, algorithm='dop853')
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


def zigzag(Z=92, E=2e18, beta=[0,1,0], trips=10, B_override=[0,0,1e-4]):
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
    interpolate_B = True
    for _ in range(2 * trips):
        print('running ' + str(_+1) + ' of ' + str(2*trips), end=': ')
        start = time.time()
        exx, eyy, epos, eb, dxx, dyy, dpos, db = run(Z=Z, E=E, betas=betas, B_override=B_override, plot=False, positions=positions, zigzag=True, interpolate_B=interpolate_B)
        print(time.time() - start)
        ex.append(exx)
        ey.append(eyy)
        dx.append(dxx)
        dy.append(dyy)
        Z = -1 * Z
        e_beta = -1. * eb
        d_beta = -1. * db
        positions = np.concatenate([epos, dpos])
        betas = np.concatenate([e_beta, d_beta])
        #interpolate_B = not interpolate_B

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
        peek = os.path.split(file)[1].split('_')[0]
        offset = 0
        if (not peek.isdigit()):
            offset = 1
        theta = float(os.path.split(file)[1].split('_')[offset])
        phi = float(os.path.split(file)[1].split('_')[offset + 1])
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
            ax.plot3D(x, y, z, color=tuple(np.random.random(3)))
            #ax.plot3D(x[:2], y[:2], z[:2], color=tuple(np.random.random(3)))
            #ax.plot(y, z, zdir='x', color=tuple(np.random.random(3)))

    gz.earth.Earth().draw(ax)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    #delta=.001
    #ax.set_xlim(1-delta, 1+delta)
    #ax.set_ylim(-delta, delta)
    #ax.set_zlim(-delta, delta)
    #ax.view_init(azim=phi-90, elev=0)
    
    
def accuracy(Z=92, E=2e18, beta=[0,1,0], divisor1=1e2, divisor2=1e3, divisor3=1e4):    
    position = gz.coordinates.Cartesian.earth
    beta = np.asarray(beta)
    beta = beta / np.sqrt(np.dot(beta,beta))
    
    out1 = gz.path.Outgoing(position, beta, Z, E, zigzag=True, save=False, max_step=.1)
    out2 = gz.path.Outgoing(position, beta, Z, E, zigzag=True, save=False, max_step=.1)
    out3 = gz.path.Outgoing(position, beta, Z, E, zigzag=True, save=False, max_step=.1)
    
    gz.path.Path.DOP853_DIVISOR = divisor1
    start1 = time.time()
    out1.propagate(B_override=[0,0,1e-4], algorithm='dop853')
    elapsed1 = time.time() - start1    
    
    gz.path.Path.DOP853_DIVISOR = divisor2
    start2 = time.time()
    out2.propagate(B_override=[0,0,1e-4], algorithm='dop853')
    elapsed2 = time.time() - start2
    
    gz.path.Path.DOP853_DIVISOR = divisor3
    start3 = time.time()
    out3.propagate(B_override=[0,0,1e-4], algorithm='dop853')
    elapsed3 = time.time() - start3

    x1 = [t[0] for t in out1.telemetry]
    y1 = [t[1] for t in out1.telemetry]
    
    x2 = [t[0] for t in out2.telemetry]
    y2 = [t[1] for t in out2.telemetry]
    
    x3 = [t[0] for t in out3.telemetry]
    y3 = [t[1] for t in out3.telemetry]

    if (plot):    
        plt.figure(figsize=[8,8])
        plt.plot(x1, y1, 'g')
        plt.plot(x2, y2, 'm')
        plt.plot(x3, y3, 'k')
        plt.show()
        
        print('Out 1 elapsed: ' + str(elapsed1))
        print('Out 2 elapsed: ' + str(elapsed2))
        print('Out 3 elapsed: ' + str(elapsed3))
        print()
        
        print('Out 1 distance: ' + str(out1.distance))
        print('Out 2 distance: ' + str(out2.distance))
        print('Out 3 distance: ' + str(out3.distance))
        print()
        
        difference = out1.position - out2.position
        difference = np.sqrt( np.dot(difference, difference) )
        print('Difference 1-2 [AU]:     ' + str(difference))
        difference *= gz.units.Change.AU_to_meter
        print('Difference 1-2 [meters]: ' + str(difference))
        print()
        
        difference = out1.position - out3.position
        difference = np.sqrt( np.dot(difference, difference) )
        print('Difference 1-3 [AU]:     ' + str(difference))
        difference *= gz.units.Change.AU_to_meter
        print('Difference 1-3 [meters]: ' + str(difference))
        print()
        
        print('out1 steps: ' + str(len(out1.telemetry)))
        print('out2 steps: ' + str(len(out2.telemetry)))
        print('out3 steps: ' + str(len(out3.telemetry)))
    