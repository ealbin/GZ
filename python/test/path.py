#!/user/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time

import gz

# a worst case, Z = 2, E = 200e18, B = 1e-12
def run(Z=92, E=2e18, B=1e-6, algorithm='dop853'):
    start = time.time()
    gyro_radius  = 1. / gz.units.SI.lightspeed * E / Z / B
    gyro_radius *= gz.units.Change.meter_to_AU
    print('Gyro-radius: ' + str(np.round(gyro_radius)) + ' [AU]')
    
    if (algorithm == 'euler'):
        step_divisor = gz.path.Path.EULER_DIVISOR
    elif (algorithm == 'dop853'):
        step_divisor = gz.path.Path.DOP853_DIVISOR
    
    position = np.asarray([gyro_radius, 0, 0])
    beta = np.asarray([0, 1., 0])
    Bfield = np.asarray([0, 0, -B])
    
    x = [position[0]]
    y = [position[1]]
    
    path = gz.path.Path(position, beta, Z, E, max_step=gyro_radius / step_divisor)
    
    r_record = [path.dist_sun]
    theta = 0.
    r2 = path.dist_sun
    closest_theta = None
    closest = None
    while (theta <= 10 * 2 * 3.2):
        path.propagate(B_override=Bfield, algorithm=algorithm)
        r_record.append(path.dist_sun)
        r1 = path.dist_sun
        step = path.step
        theta += np.arccos( (r1**2 + r2**2 - step**2) / (2. * r1 * r2) )
        r2 = r1
        x.append(path.position[0])
        y.append(path.position[1])
        
        if (theta > 2. * 3.):
            if (closest_theta is None):
                closest_theta = theta
                closest = r1
            elif (np.abs(theta - 2.*np.pi) < np.abs(closest_theta - 2.*np.pi)):
                closest_theta = theta
                closest = r1
            
    plt.figure(figsize=[8,8])
    plt.plot(x, y)
    plt.xlim(-1.1*gyro_radius, 1.1*gyro_radius)
    plt.ylim(-1.1*gyro_radius, 1.1*gyro_radius)
    plt.grid()
    plt.show()
    
    print('Finish distance - gyro-radius [AU]: ' + str(np.abs(closest - gyro_radius)))
    print('theta [deg]: ' + str(180./np.pi * closest_theta))
    print("done in [sec]: " + str(time.time() - start))
    print('last stepsize: ' + str(path.step))
    
    r_record = np.asarray(r_record)
    mean = np.mean(r_record)
    stddev = np.std(r_record)
    print('Mean radius [AU]: ' + str(mean))
    print('StdDev [AU]: ' + str(stddev))
    print('StdDev [m]:  ' + str(stddev * gz.units.Change.AU_to_meter))