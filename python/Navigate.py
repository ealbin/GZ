
import numpy as np
import LorentzForce

def newHeading(position, beta, ratio, distance):
    beta = np.array(beta)
    new_beta = beta + LorentzForce.dBeta_dS(position, beta, ratio) * distance
    new_beta /= np.sqrt( np.dot(new_beta, new_beta) )
    return new_beta

def stepPosition(position, beta, ratio, distance):
    position = np.array(position)
    new_position = position + newHeading(position, beta, ratio, distance) * distance
    return new_position
