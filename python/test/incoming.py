#!/user/bin/env python3

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time

import gz

def run(Z=92, E=2e18, beta=[0, 1, 0], B_override=[0,0,1e-4],
        position=gz.coordinates.Cartesian.earth, seed=12345, max_step=.1, 
        delta=50 * gz.units.SI.radius_earth * gz.units.Change.meter_to_AU):
    start = time.time()
    
    beta = np.asarray(beta)
    beta = beta / np.sqrt(np.dot(beta, beta))
    
    outgoing = gz.path.Outgoing(position, beta, Z, E, save=False) #, max_step=max_step)
    outgoing.propagate(B_override=B_override)
    
    origin = position
    position = outgoing.position
    beta = -1. * outgoing.beta
    A = gz.units.Nuclide.mass_number(Z)
    dist = outgoing.distance / 4.
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

    fig = plt.figure(figsize=[16,16])
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
    
    fig.show()
