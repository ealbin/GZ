#!/user/bin/env python3

import numpy as np
import gz

def run(E_cr, Z_cr, Beta_cr, E_ph, Beta_ph):

    A_cr = gz.units.Nuclide.mass_number(Z_cr)
    M_cr = A_cr * gz.units.Change.amu_to_eV
    M_ph = 0.
    
    Beta_cr = np.asarray(Beta_cr)
    Beta_ph = np.asarray(Beta_ph)
    Beta_cr = Beta_cr / np.sqrt(np.dot(Beta_cr, Beta_cr))
    Beta_ph = Beta_ph / np.sqrt(np.dot(Beta_ph, Beta_ph))
    
    P_cr = gz.relativity.momentum(E_cr, M_cr, Beta_cr)
    P_ph = gz.relativity.momentum(E_ph, M_ph, Beta_ph)
    
    E_net = E_cr + E_ph
    P_net = P_cr + P_ph
    
    ####
    
    E_d = (A_cr - 1.) / A_cr * E_net
    E_n =         1.  / A_cr * E_net
    
    M_d = (A_cr - 1.) * gz.units.Change.amu_to_eV
    M_n =         1.  * gz.units.Change.amu_to_eV
    
    P_d_mag = gz.relativity.momentum_mag(E_d, M_d)
    P_n_mag = gz.relativity.momentum_mag(E_n, M_n)
    
    numerator = np.dot(P_net, P_net) - P_d_mag**2 - P_n_mag**2
    denominator = 2. * P_d_mag * P_n_mag
    
    cosTheta = numerator / denominator
    if (cosTheta > 1. and np.isclose(cosTheta, 1.)):
        cosTheta = 1.
    return np.arccos(cosTheta) * 180. / np.pi
        