"""Compute the interaction cross section for photodissintegration
"""
import numpy as np

###########################################################################    

### parametric constants:

###########################################################################    


def lorentzBackBoost( En_eV, mass, angle_rad ):
    """Returns the original energy of a photon measured from the reference frame of the Sun.
    That is, back-boosted from the energy encountered in the frame of a traveling nucleus En_eV 
    in [electronVolts] with mass [atomic mass units].  The relative angle between
    photon and nucleus momenta is given by angle_rad in [radians].
    Head-on, angle_rad = 0; same direction, angle_rad = pi.
    ref. 1960 Geramsimova and Zatsepin
    """
    amu2mev = 931.5 # [MeV/c**2] e.g. carbon has an atomic mass of 12, or 12 * amu2mev = 11,178 MeV/c**2
    mass_eV = mass * amu2mev * 1e6 # [eV / c**2]
    
    lorentz_gamma = En_eV / mass_eV
    
    return En_eV / ( 2. * lorentz_gamma * np.cos(angle_rad / 2.)**2 )

def solarPhotonDensity( rs_AU, Es_eV ):
    """Returns the differential solar photon number density dn/dE in 
    [number / eV * cm**3] given a radial distance from the Sun, rs_AU in
    [astronomical units] and solar photon energy Es_eV in [electronVolts] 
    as measured in the reference frame of the Sun.
    Black body spectrum with T = 5770 K.
    """
    return 7.8e7 * (1./rs_AU**2) * ( Es_eV**2 / ( np.exp(Es_eV/.5) - 1. ) )

def singleNucleon( mass_number, photon_energy_eV ):
    """Returns the photodisintegration cross section for losing one nucleon by a nucleus 
    of mass_number [unit-less] (a.k.a. "A") through interaction with a photon with energy
    photon_energy [eV] in the nucleus' frame of reference.
    Reference 1999 Epele, Mollerach and Roulet.
    """
    def giantDipoleResonance( A, E_MeV ):
        """Returns cross section model for GDR interaction [cm**2].
        """
        sigma0 = 1.45e-27 * A # [cm**2], cross section scale factor
        T      = 8.           # [MeV], GDR energy bandwidth
        if A <= 4:
            epsilon0 = 0.925 * A**2.433 # [MeV], peak energy of GDR resonance
        else:
            epsilon0 = 42.65 * A**-0.21 # [MeV]
        
        numerator   = ( E_MeV * T )**2    
        denominator = ( E_MeV**2 - epsilon0**2 )**2 + ( E_MeV * T )**2
        shapefactor = numerator / float(denominator) # [unit-less] peak shape factor
        
        return sigma0 * shapefactor # [cm**2]
    
    def prePionProduction( A, E_MeV ):
        """Returns cross section model for energies between 30 and 150 MeV.
        note: quasi-deuteron or multiple nucleon ejection turns on in this regeme.
        This cross section represents single nucleon ejection only.  Proceed with caution.
        """
        low_bound = A / 8. * 1e-27 # [cm**2]
        gdr_bound = giantDipoleResonance( A, E_MeV ) # [cm**2]
        return max([ gdr_bound, low_bound ])

    def postPionProduction( A, E_MeV ):
        """Returns cross section model for energies above 150 MeV.
        note: pion production turns on around 150 MeV and nucleons are knocked out
        via photon-absorption with nearest resonance at the Delta baryon mass 1232 MeV
        (proton ~938 MeV + ~300 MeV photon).
        Multiple nucleon emission is increasing likely.. use caution with results.
        """
        S  = 0.3
        nu = 1.8
        epsilon1  = 180 # [MeV]
        epsilon_t = ( E_MeV - 150. ) / epsilon1
        
        piece_1 = A / 8.
        piece_2 = A * S * epsilon_t * np.exp( (1 - epsilon_t**nu) / nu )
        return (piece_1 + piece_2) * 1e-27 # [cm**2]
                
    
    photon_energy_MeV = photon_energy_eV / 1e6
    
    if photon_energy_MeV <= 30.:  # i.e. 30 [MeV]
        return giantDipoleResonance( mass_number, photon_energy_MeV )
    elif photon_energy_MeV <= 150.: # i.e. 150 [MeV]
        return prePionProduction( mass_number, photon_energy_MeV )
    else:
        return postPionProduction( mass_number, photon_energy_MeV )
        
        