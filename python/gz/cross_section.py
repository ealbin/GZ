#!/usr/bin/env python3

"""Compute the interaction cross section for photodissintegration
"""

__project__     = 'GZ Paper'
__version__     = 'v1.0'
__objective__   = 'Phenominology'
__institution__ = 'University of California, Irvine'
__department__  = 'Physics and Astronomy'
__author__      = 'Eric Albin'
__email__       = 'Eric.K.Albin@gmail.com'
__updated__     = '13 May 2019'


import numpy as np

from . import units


class Photodissociation:
    
    # TODO: use nuclear data instead of a model

    def singleNucleon(proton_number, photon_energy_eV, mass_number=None):
        """Returns the photodisintegration cross section [cm**2] for losing one nucleon by a nucleus 
        of mass_number [unit-less] (a.k.a. "A") through interaction with a photon with energy
        photon_energy [eV] in the nucleus' frame of reference.
        Reference 1999 Epele, Mollerach and Roulet.
        If mass_number is None, uses average mass number from units module.
        """
        
        if (mass_number == None):
            mass_number = units.Nuclide.mass_number(proton_number)
        
        
        def giantDipoleResonance(A, E_MeV):
            """Returns cross section model for GDR interaction [cm**2].
            """
            sigma0 = 1.45e-27 * A # [cm**2], cross section scale factor
            T      = 8.           # [MeV], GDR energy bandwidth
            if (A <= 4):
                epsilon0 = 0.925 * A**2.433 # [MeV], peak energy of GDR resonance
            else:
                epsilon0 = 42.65 * A**-0.21 # [MeV]
            
            numerator   = (E_MeV * T)**2    
            denominator = (E_MeV**2 - epsilon0**2)**2 + (E_MeV * T)**2
            shapefactor = numerator / float(denominator) # [unit-less] peak shape factor
            
            return sigma0 * shapefactor # [cm**2]
        
        
        def prePionProduction(A, E_MeV):
            """Returns cross section model for energies between 30 and 150 MeV.
            note: quasi-deuteron or multiple nucleon ejection turns on in this regeme.
            This cross section represents single nucleon ejection only.  Proceed with caution.
            """
            low_bound = A / 8. * 1e-27 # [cm**2]
            gdr_bound = giantDipoleResonance(A, E_MeV) # [cm**2]
            return max([gdr_bound, low_bound])
    
    
        def postPionProduction(A, E_MeV):
            """Returns cross section model for energies above 150 MeV.
            note: pion production turns on around 150 MeV and nucleons are knocked out
            via photon-absorption with nearest resonance at the Delta baryon mass 1232 MeV
            (proton ~938 MeV + ~300 MeV photon).
            Multiple nucleon emission is increasing likely.. use caution with results.
            """
            S  = 0.3
            nu = 1.8
            epsilon1  = 180 # [MeV]
            epsilon_t = (E_MeV - 150.) / epsilon1
            
            piece_1 = A / 8.
            piece_2 = A * S * epsilon_t * np.exp( (1 - epsilon_t**nu) / nu )
            return (piece_1 + piece_2) * 1e-27 # [cm**2]
                    
        
        photon_energy_MeV = photon_energy_eV / 1e6
        
        if (photon_energy_MeV <= 30.):  # i.e. 30 [MeV]
            return giantDipoleResonance(mass_number, photon_energy_MeV)
        elif (photon_energy_MeV <= 150.): # i.e. 150 [MeV]
            return prePionProduction(mass_number, photon_energy_MeV)
        else:
            return postPionProduction(mass_number, photon_energy_MeV)
        
