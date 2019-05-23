#!/user/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time

import gz

def run(Z=92, E=2e18, beta=[0, 1, 0], B_override=[0,0,1e-4],
        position=gz.coordinates.Cartesian.earth, seed=12345, max_step=.01, 
        delta=50 * gz.units.SI.radius_earth * gz.units.Change.meter_to_AU):
    start = time.time()
    
    beta = np.asarray(beta)
    beta = beta / np.sqrt(np.dot(beta, beta))
    
    outgoing = gz.path.Outgoing(position, beta, Z, E, save=False, max_step=max_step)
    outgoing.propagate(B_override=B_override)
    
    origin = position
    position = outgoing.position
    beta = -1. * outgoing.beta
    A = gz.units.Nuclide.mass_number(Z)

    dists = []
    probs = []
    max_dist = outgoing.telemetry[-1][6]
    length = len(outgoing.telemetry)
    for _ in range(length):
        t = outgoing.telemetry[length - _ - 1]
        pos = t[:3]
        bet = -1. * t[3:6]
        dis = max_dist - t[6]
        atten = gz.probability.Solar.attenuation(pos, bet, Z, E, mass_number=A)
        probs.append(gz.probability.oneOrMore(atten, dis))
        dists.append(dis)    
    dist = gz.probability.random(dists, probs, 1, seed=seed)   
    
    incoming = gz.path.Incoming(origin, position, beta, Z, A, E, dist, max_step=max_step, save=False)
    incoming.propagate(B_override=B_override, seed=seed)
    
    elapsed_time = time.time() - start

    out_x = [t[0] for t in outgoing.telemetry]
    out_y = [t[1] for t in outgoing.telemetry]
    out_z = [t[2] for t in outgoing.telemetry]
    
    in_x = [t[0] for t in incoming.telemetry]
    in_y = [t[1] for t in incoming.telemetry]
    in_z = [t[2] for t in incoming.telemetry]
    
    p_x = [t[0] for t in incoming.p_path.telemetry]
    p_y = [t[1] for t in incoming.p_path.telemetry]
    p_z = [t[2] for t in incoming.p_path.telemetry]
    n_x = [t[0] for t in incoming.n_path.telemetry]
    n_y = [t[1] for t in incoming.n_path.telemetry]
    n_z = [t[2] for t in incoming.n_path.telemetry]

    dp_x = [t[0] for t in incoming.dp_path.telemetry]
    dp_y = [t[1] for t in incoming.dp_path.telemetry]
    dp_z = [t[2] for t in incoming.dp_path.telemetry]
    dn_x = [t[0] for t in incoming.dn_path.telemetry]
    dn_y = [t[1] for t in incoming.dn_path.telemetry]
    dn_z = [t[2] for t in incoming.dn_path.telemetry]

    def dist_earth(x, y, z):
        pos = np.asarray([x, y, z]).T - gz.coordinates.Cartesian.earth
        dist = np.zeros(len(pos))
        for _ in range(len(pos)):
            dist[_] = np.sqrt(np.dot(pos[_], pos[_]))
        dist = dist / ( gz.units.SI.radius_earth * gz.units.Change.meter_to_AU )
        return dist

    o_dist = dist_earth(out_x, out_y, out_z)
    i_dist = dist_earth(in_x, in_y, in_z)
    p_dist = dist_earth(p_x, p_y, p_z)
    n_dist = dist_earth(n_x, n_y, n_z)
    dp_dist = dist_earth(dp_x, dp_y, dp_z)
    dn_dist = dist_earth(dn_x, dn_y, dn_z)

    fig1 = plt.figure(figsize=[16,16])
    x = [y for y in range(len(o_dist))]
    start = x[-1]
    plt.plot(x, o_dist, 'gray')
    x = [y for y in range(start, start + len(i_dist))]
    plt.plot(x, i_dist, 'k')
    start = x[-1]
    x = [y for y in range(start, start + len(p_dist))]
    plt.plot(x, p_dist, 'b')
    x = [y for y in range(start, start + len(n_dist))]
    plt.plot(x, n_dist, 'r')
    x = [y for y in range(start, start + len(dp_dist))]
    plt.plot(x, dp_dist, 'm')
    x = [y for y in range(start, start + len(dn_dist))]
    plt.plot(x, dn_dist, 'y')
    plt.yscale('log')
    plt.show()
    plt.ylim(1e-1, 10 * max(np.concatenate([o_dist, i_dist, p_dist, n_dist, dp_dist, dn_dist])))
    
    fig2 = plt.figure(figsize=[16,16])
    ax = plt.axes(projection='3d')
    ax.plot3D(out_x, out_y, out_z, 'gray')
    ax.plot3D(in_x, in_y, in_z, 'k')
    ax.plot3D(p_x, p_y, p_z, 'b')
    ax.plot3D(dp_x, dp_y, dp_z, 'm')
    ax.plot3D(n_x, n_y, n_z, 'r')
    ax.plot3D(dn_x, dn_y, dn_z, 'y')
    
    # draw earth sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v) * gz.units.SI.radius_earth * gz.units.Change.meter_to_AU
    y = np.sin(u)*np.sin(v) * gz.units.SI.radius_earth * gz.units.Change.meter_to_AU
    z = np.cos(v) * gz.units.SI.radius_earth * gz.units.Change.meter_to_AU
    x += gz.coordinates.Cartesian.earth[0]
    y += gz.coordinates.Cartesian.earth[1]
    z += gz.coordinates.Cartesian.earth[2]
    ax.plot_wireframe(x, y, z, color='gray')
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.set_xlim(1-delta, 1+delta)
    ax.set_ylim(-delta, delta)
    ax.set_zlim(-delta, delta)    
    ax.view_init(azim=-45, elev=-10)    
    
    print('Elapsed: ' + str(elapsed_time))
    print('gray=outgoing, black=incoming')
    print('blue=proton, magenta=daughter')
    print('red=neutron, yellow=daughter')
    print()
    print('Final distance proton:   ' + str(p_dist[-1]))
    print('Final distance neutron:  ' + str(n_dist[-1]))
    print('Final distance dproton:  ' + str(dp_dist[-1]))
    print('Final distance dneutron: ' + str(dn_dist[-1]))
