import time
import numpy as np
import LorentzForce

# all distance parameters in AU

AU2M = 149597870700. # [AU] * au2m = [m]

def newHeading(position, beta, ratio, stepsize):
    beta = np.array(beta)
    new_beta = beta + LorentzForce.dBeta_dS(position, beta, ratio) * stepsize * AU2M
    new_beta /= np.sqrt( np.dot(new_beta, new_beta) )
    return new_beta

def stepPosition(position, beta, ratio, stepsize):
    position = np.array(position)
    new_position = position + newHeading(position, beta, ratio, stepsize) * stepsize
    return new_position

def runCourse(position, beta, ratio, stepsize, AU_limit=5, filename='telemetry.txt'):
    start    = time.time()
    count    = 0
    history  = []
    position = np.array(position)
    while np.dot(position, position) < AU_limit**2:
        history.append(position)
        position = stepPosition(position, beta, ratio, stepsize)
        count += 1
    history.append(position)
    end = time.time()
    elapsed = end - start
    print 'total elapsed time: {}, ave time/step: {}, steps: {}'.format(elapsed, elapsed/float(count), count)
    print 'final coordinates: {}'.format(position)
    with open(filename, 'w') as file:
        for point in history:
            file.write('{}\n'.format(point))
