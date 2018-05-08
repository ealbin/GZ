# returns d(beta)/ds from the magnetic field interaction

import numpy as np
import MagneticModel 

c = 3.e8 # m/s speed of light

# returns new beta after interacting with the magnetic field
def dBeta_dS( position, beta, ratio ):
    position = np.array(position)
    beta     = np.array(beta)

    B = MagneticModel.Bfield(position)
    
    return ratio * np.cross( c*beta, B )
