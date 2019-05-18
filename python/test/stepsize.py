#!/user/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import time

def run(divisor=10., radius=10.):
    start = time.time()
    
    step = radius / divisor
    position = np.asarray([radius, 0, 0])
    beta = np.asarray([0, 1., 0])
    theta = 0.
    
    distance = []
    x = []
    y = []
    distance.append(np.sqrt(np.dot(position, position)))
    r2 = distance[0]
    x.append(position[0])
    y.append(position[1])
    while (theta <= 2 * 3.2):
    
        beta = beta / np.sqrt(np.dot(beta,beta))    
        unit_v = -np.cross(beta, np.asarray([0,0,1]))
        dbeta_ds = 1. / radius * unit_v
        
        beta += dbeta_ds * step
        position += beta * step
        r1 = np.sqrt(np.dot(position, position))
        if (len(distance) < 3):
            distance.append(r1)
        #x.append(position[0])
        #y.append(position[1])
        theta += np.arccos( (r1**2 + r2**2 - step**2) / (2. * r1 * r2) )
        r2 = r1
    
    distance = np.asarray(distance)
    x = np.asarray(x)
    y = np.asarray(y)
    
    plt.figure(figsize=[8,8])
    plt.plot(x,y)
    plt.xlim(-1.2*radius, 1.2*radius)
    plt.ylim(-1.2*radius, 1.2*radius)
    plt.grid()
    plt.show()
    
    plt.figure(figsize=[8,8])
    plt.plot(distance)
    plt.ylim(radius - step, radius + step)
    plt.grid()
    plt.show()
    
    print(np.abs(distance[1] - radius))
    print('elapsed time [sec]: ' + str(time.time() - start))
    return x, y, distance