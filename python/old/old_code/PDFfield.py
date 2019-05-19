#!/usr/bin/env python

"""Precompute the probability density field for decay, store to disk 
and use it as an interpolated look-up table to profoundly accelerate
numeric integration.
"""

import numpy as np
import os
import sys
import tarfile
import time

from scipy import interpolate
from scipy import integrate
from scipy import optimize

import Constants
import PDFmodel

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Development"


# global field values in memory
#-------------------------------
__spacelimit = None
__resolution = None
__x   = None
__y   = None
__z   = None
__PDF = None
__InterpolatePDF = None

def precompute(mass_number, energy_eV, spacelimit=6, resolution=60, autoload=True, directory='tables'):
    """Returns probability density field x, y, z, PDF meshes by
    by disk-read or re-generation.  Probility density is in [probability / meters].
    
    mass_number: atmoic mass "A" of isotope (number of protons plus neutrons)
    
    energy_eV  : energy of isotope [electronVolts]
    
    spacelimit : radial reach (r) of the space volume 
                (x, y, z) === (-r to r) by (-r to r) by (-r to r) [astronomical units]

    resolution : the number of samples taken between (-r to r) along each dimension.
                 In addition, there are another resolution's-worth of samples added to
                 that set between (-r/10 to r/10) to resolve near the Sun better.
                 resolution = 60 takes around (TODO) hours to regenerate.
    
    autoload   : if True, look FIRST to disk for existing table.
                    if (no preexisting) or (has different spacelimit or resolution):
                        regenerate from scratch and overwrite existing.
                 if False, force regenerate from scratch and overwrite existing.
    
    directory  : subdirectory with magnetic field text file
    
    returns dictionary { 'x', 'y', 'z', 'PDF' }
                 x, y, z and PDF have shape (<=2*resolution,)
                 The <=2 is because some points are common to both (-r to r) and 
                 (-r/10 to r/10), thus the shape is between (1 to 2)*resolution.

    """
    # check if already loaded in memory, return and exit if so
    #----------------------------------------------------------
    global __spacelimit, __resolution
    global __x, __y, __z, __PDF
    
    if ( (autoload) and (__spacelimit == spacelimit) and (__resolution == resolution) and
         (__x is not None)  and (__y is not None)  and (__z is not None) and (__PDF is not None) ):
        return {'x':__x, 'y':__y, 'z':__z, 'PDF':__PDF}    

    
    # configure x,y,z and pdf
    #------------------------
    spacelimit = int( spacelimit ) # [astronomical units] (integer for easy file read)
    resolution = int( resolution ) # N divisions (integer for easy file read)

    x = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    y = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    z = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    ## add extra points around the sun:
    x = np.union1d(x, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    y = np.union1d(y, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    z = np.union1d(z, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    ## add extra points in the sun-earth corridor:
    x = np.union1d(x, np.linspace(    0,    1, resolution) )
    y = np.union1d(y, np.linspace(-0.15, 0.15, resolution) )
    z = np.union1d(z, np.linspace(-0.15, 0.15, resolution) )
    
    spacelimit = np.array([spacelimit])
    resolution = np.array([resolution])

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    shape = np.array(X.shape)
    size = X.flatten().size
    
    pdf = np.zeros(size) # [probability / meter]

    # load from file or regenerate
    #------------------------------
    regen = False
    if not autoload:
        regen = True

    text_sep = ', '
    base_dir  = os.path.dirname( os.path.abspath( __file__ ) ) 
    directory = 'tables'
    p_fname   = 'A_{0}_E_{1}.pdf'.format(mass_number, energy_eV)
    p_fnameZip= p_fname + '.tar.gz' 
    p_path    = os.path.abspath( os.path.join( base_dir, directory, p_fname ) )
    p_exists  = os.path.isfile(p_path)
    if (not p_exists):
        zip_path = os.path.abspath( os.path.join( base_dir, directory, p_fnameZip ) )
        if os.path.isfile(zip_path):
            with tarfile.open(zip_path, 'r:gz') as tf:
                tf.extractall( os.path.abspath( os.path.join( base_dir, directory) ) )
                return precompute(spacelimit=spacelimit, resolution=resolution, autoload=autoload, directory=directory, p_fname=p_fname)
        else:
            regen = True

    # load from file if the file is good
    if (not regen):
        with open(p_path) as p_f:
            header = p_f.readline()
            f_spacelimit = int( p_f.readline() )
            f_resolution = int( p_f.readline() )
            f_shape      = np.fromstring( p_f.readline(), sep=text_sep )                            
            if ( f_spacelimit == spacelimit and 
                 f_resolution == resolution and 
                 f_shape == shape).all():
                x   = np.fromstring( p_f.readline(), sep=text_sep)
                y   = np.fromstring( p_f.readline(), sep=text_sep)
                z   = np.fromstring( p_f.readline(), sep=text_sep)
                pdf = np.fromstring( p_f.readline(), sep=text_sep)
            else:
                regen = True

    # regenerate and overwrite
    if regen:
        i_max  = X.flatten().size
        target = 0.
        start  = time.time()
        for i, (ix, iy, iz) in enumerate( zip( X.flatten(), Y.flatten(), Z.flatten() ) ):
            pdf[i] = PDFmodel.pdf( np.asarray([ix, iy, iz]), mass_number, energy_eV ) # [probability / meter]

            # progress report for long regenerations
            if ( i / float(i_max) ) >= ( target / 100. ):
                print('\r                                                           \r',)
                print('  progress: {:.1f}%   elapsed: {:.2f} [sec]'.format(target, time.time() - start ),)
                sys.stdout.flush()
                target += .1
        print
        
        with open(p_path, 'w') as p_f:
            header = ( 'rows: 0:this header, 1:spacelimit [AU], 2:resolution, 3:shape, '
                       '4:x [AU], 5:y [AU], 6:z [AU], 7:PDF [probability/meter]' )

            # header and parameters
            p_f.write(header + '\n')
            spacelimit.tofile(p_f, sep=text_sep)
            p_f.write('\n')
            resolution.tofile(p_f, sep=text_sep)
            p_f.write('\n')
            shape.tofile(p_f, sep=text_sep)
            p_f.write('\n')

            # x, y, z
            x.tofile(p_f, sep=text_sep)
            p_f.write('\n')
            y.tofile(p_f, sep=text_sep)
            p_f.write('\n')
            z.tofile(p_f, sep=text_sep)
            p_f.write('\n')

            # pdf
            pdf.tofile(p_f, sep=text_sep)
            p_f.write('\n')
            
            #### OPTIONAL TODO:
            # make tar.gz file
        
    # load into memory
    #------------------
    PDF = pdf.reshape(shape)

    __spacelimit = spacelimit
    __resolution = resolution
    __x   = x
    __y   = y
    __z   = z
    __PDF = PDF
    
    return {'x':__x, 'y':__y, 'z':__z, 'PDF':__PDF}    


def cartesianPDF( cartesian_pos, mass_number, energy_eV, close2sun=0.01 ):
    """Returns cartesian [probability / meter] value (PDF) for cartesian_pos = (x, y, z).
    If position is within close2sun radius [AU], do not interpolate, return exact (slower).
    For spacelimit==6 and resolution==60, interpolation is acceptable up to close2sun==0.01.
    """
    cartesian_pos = np.array(cartesian_pos)
    if np.sqrt(np.dot(cartesian_pos, cartesian_pos)) < close2sun:
        return PDFmodel.pdf(cartesian_pos, mass_number, energy_eV)
    
    global __InterpolatePDF
    if (__InterpolatePDF is not None):
        try:
            PDF = __InterpolatePDF(cartesian_pos)
        except:
            # out of bounds, compute explicitly
            PDF = PDFmodel.pdf(cartesian_pos, mass_number, energy_eV)
        return PDF
    else:
        meshes = precompute(mass_number, energy_eV)
        x   = meshes['x']
        y   = meshes['y']
        z   = meshes['z']
        PDF = meshes['PDF']
        
        __InterpolatePDF = interpolate.RegularGridInterpolator((x,y,z), PDF, bounds_error=True) # [Probability / meter]
        
        return cartesianPDF(cartesian_pos, mass_number, energy_eV)

    
def integratePath( theta_x, phi_x, r, mass_number, energy_eV ):
    """ Returns the integrated probability for a path specified by theta_x and phi_x, 
    decaying a distance r away from Earth.  The nuclei has mass_number and energy_eV.
    Returns dictionary of 'probability' (integrated) and 'path' (path-distance, probability).
    """
    # integration parameters
    stepsize_min = Constants.meters2AU(1e8 ) # [AU]
    stepsize_max = Constants.meters2AU(1e11) # [AU]
    integration_limit = 15 # [AU]
    dist_min = 1. # switch to stepsize_min [AU]
    dist_max = 8. # switch to stepsize_max [AU]
    
    tx_rad    = theta_x * np.pi / 180. # [radians]
    px_rad    = phi_x   * np.pi / 180. # [radians]
    heading_x = np.cos(tx_rad)
    heading_y = np.sin(tx_rad) * np.cos(px_rad)
    heading_z = np.sin(tx_rad) * np.sin(px_rad)
    heading   = np.array([ heading_x, heading_y, heading_z ])
    
    earth_pos = np.array([1, 0, 0]) # [AU]
    pos = earth_pos + r * heading   # [AU]
    
    # for converting probability / meter into probability / AU
    m_per_AU = Constants.AU2meters(1.) # [meters / AU]
    mass_number = int(mass_number)
    
    dists = [ 0 ]
    pdfs  = [ cartesianPDF(pos, mass_number, energy_eV) * m_per_AU ]
    while np.linalg.norm(pos) < integration_limit:
        solar_dist = np.linalg.norm(pos)
        if solar_dist < dist_min: 
            stepsize = stepsize_min
        elif solar_dist > dist_max:
            stepsize = stepsize_max
        else:
            stepsize = (solar_dist - dist_min) * (stepsize_max - stepsize_min) / (dist_max - dist_min) + stepsize_min
        pos = pos + stepsize * heading
        dists.append( dists[-1] + stepsize )
        pdfs.append( cartesianPDF(pos, mass_number, energy_eV) * m_per_AU )
    
    probability = integrate.simps(pdfs, x=dists, even='first')
    return {'probability':probability, 'path':np.array([dists, pdfs])}

