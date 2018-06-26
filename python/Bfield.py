#!/usr/bin/env python

"""Precompute the total magnetic field, store to disk 
and use it as an interpolated look-up table to profoundly accelerate
numeric integration.
"""

import numpy as np
import os
import sys
import tarfile
import time

from scipy import interpolate
from scipy import optimize

import SolarMagneticModel
import Transform

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Production"

# global field values in memory
#-------------------------------
__spacelimit = None
__resolution = None
__x  = None
__y  = None
__z  = None
__BX = None
__BY = None
__BZ = None
__InterpolateBx = None
__InterpolateBy = None
__InterpolateBz = None

def precompute(spacelimit=6, resolution=60, autoload=True, directory='tables', b_fname='cartesianBfield.Tesla'):
    """Returns total magnetic field x, y, z, BX, BY, BZ meshes by
    by disk-read or re-generation.  Field density in Teslas.
    
    spacelimit : radial reach (r) of the space volume 
                (x, y, z) === (-r to r) by (-r to r) by (-r to r) [astronomical units]

    resolution : the number of samples taken between (-r to r) along each dimension.
                 In addition, there are another resolution's-worth of samples added to
                 that set between (-r/10 to r/10) to resolve near the Sun better.
                 resolution = 60 takes around 5 hours to regenerate.
    
    autoload   : if True, look FIRST to disk for existing table.
                    if (no preexisting) or (has different spacelimit or resolution):
                        regenerate from scratch and overwrite existing.
                 if False, force regenerate from scratch and overwrite existing.
    
    directory  : subdirectory with magnetic field text file
    
    b_fname    : filename for magnetic field text file
    
    returns dictionary { 'x', 'y', 'z', 'BX', 'BY', 'BZ' }
                 x, y, z have shape (<=2*resolution,)
                 BX, BY, BZ have shape (<=2*resolution, <=2*resolution, <=2*resolution)
                 The <=2 is because some points are common to both (-r to r) and 
                 (-r/10 to r/10), thus the shape is between (1 to 2)*resolution.

    """
    # check if already loaded in memory, return and exit if so
    #----------------------------------------------------------
    global __spacelimit, __resolution
    global __x, __y, __z, __BX, __BY, __BZ
    
    if ( (__spacelimit == spacelimit) and (__resolution == resolution) and
         (__x is not None)  and (__y is not None)  and (__z is not None) and
         (__BX is not None) and (__BY is not None) and (__BZ is not None) ):
        return {'x':__x, 'y':__y, 'z':__z, 'BX':__BX, 'BY':__BY, 'BZ':__BZ}    

    
    # configure x,y,z and bx,by,bz
    #------------------------------
    spacelimit = int( spacelimit ) # [astronomical units] (integer for easy file read)
    resolution = int( resolution ) # N divisions (integer for easy file read)

    x = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    y = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    z = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    ## add extra points around the sun:
    x = np.union1d(x, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    y = np.union1d(y, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    z = np.union1d(z, np.linspace(-spacelimit/10., spacelimit/10., resolution) )
    
    spacelimit = np.array([spacelimit])
    resolution = np.array([resolution])

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    shape = np.array(X.shape)
    size = X.flatten().size
    
    bx = np.zeros(size) # [Tesla]
    by = np.zeros(size) # [Tesla]
    bz = np.zeros(size) # [Tesla]    

    # load from file or regenerate
    #------------------------------
    regen = False
    if not autoload:
        regen = True

    text_sep = ', '
    base_dir  = os.path.dirname( os.path.abspath( SolarMagneticModel.__file__ ) )
    directory = 'tables'
    b_fname   = 'cartesianBfield.Tesla'
    b_fnameZip= b_fname + '.tar.gz' 
    b_path    = os.path.abspath( os.path.join( base_dir, directory, b_fname ) )
    b_exists  = os.path.isfile(b_path)
    if (not b_exists):
        zip_path = os.path.abspath( os.path.join( base_dir, directory, b_fnameZip ) )
        if os.path.isfile(zip_path):
            with tarfile.open(zip_path, 'r:gz') as tf:
                tf.extractall( os.path.abspath( os.path.join( base_dir, directory) ) )
                return precompute(spacelimit=spacelimit, resolution=resolution, autoload=autoload, directory=directory, b_fname=b_fname)
        else:
            regen = True
    
    # load from file if the file is good
    if (not regen):
        with open(b_path) as b_f:
            header = b_f.readline()
            f_spacelimit = int( b_f.readline() )
            f_resolution = int( b_f.readline() )
            f_shape      = np.fromstring( b_f.readline(), sep=text_sep )                            
            if ( f_spacelimit == spacelimit and 
                 f_resolution == resolution and 
                 f_shape == shape).all():
                x  = np.fromstring( b_f.readline(), sep=text_sep)
                y  = np.fromstring( b_f.readline(), sep=text_sep)
                z  = np.fromstring( b_f.readline(), sep=text_sep)
                bx = np.fromstring( b_f.readline(), sep=text_sep)
                by = np.fromstring( b_f.readline(), sep=text_sep)
                bz = np.fromstring( b_f.readline(), sep=text_sep)
            else:
                regen = True

    # regenerate and overwrite
    if regen:
        i_max  = X.flatten().size
        target = 0.
        start  = time.time()
        for i, (ix, iy, iz) in enumerate( zip( X.flatten(), Y.flatten(), Z.flatten() ) ):
            b_solar = SolarMagneticModel.sumBfieldTesla( np.array([ ix, iy, iz ]) )
            bx[i] = b_solar[0] # [Tesla]
            by[i] = b_solar[1] # [Tesla]
            bz[i] = b_solar[2] # [Tesla]
            # progress report for long regenerations
            if ( i / float(i_max) ) >= ( target / 100. ):
                print '\r                                                           \r',
                print '  progress: {:.1f}%   elapsed: {:.2f} [sec]'.format(target, time.time() - start ),
                sys.stdout.flush()
                target += .1
        print
        
        with open(b_path, 'w') as b_f:
            header = ( 'rows: 0:this header, 1:spacelimit [AU], 2:resolution, 3:shape, '
                       '4:x [AU], 5:y [AU], 6:z [AU], 7:BX [T], 8:BY [T], 9:BZ [T]' )

            # header and parameters
            b_f.write(header + '\n')
            spacelimit.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            resolution.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            shape.tofile(b_f, sep=text_sep)
            b_f.write('\n')

            # x, y, z
            x.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            y.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            z.tofile(b_f, sep=text_sep)
            b_f.write('\n')

            # bx, by, bz
            bx.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            by.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            bz.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            #### TODO:
            # make tar.gz file
        
    # load into memory
    #------------------
    BX = bx.reshape(shape)
    BY = by.reshape(shape)
    BZ = bz.reshape(shape)

    __spacelimit = spacelimit
    __resolution = resolution
    __x  = x
    __y  = y
    __z  = z
    __BX = BX
    __BY = BY
    __BZ = BZ
    
    return {'x':__x, 'y':__y, 'z':__z, 'BX':__BX, 'BY':__BY, 'BZ':__BZ}    


def cartesianTesla( cartesian_pos, close2sun=0.01 ):
    """Returns cartesian [Tesla] values (Bx, By, Bz) for cartesian_pos = (x, y, z).
    If position is within close2sun radius [AU], do not interpolate, return exact (slow).
    For spacelimit==6 and resolution==60, interpolation is acceptable up to close2sun==0.01.
    """
    cartesian_pos = np.array(cartesian_pos)
    if np.sqrt(np.dot(cartesian_pos, cartesian_pos)) < close2sun:
        return SolarMagneticModel.sumBfieldTesla(cartesian_pos)
    
    global __InterpolateBx, __InterpolateBy, __InterpolateBz
    if ( (__InterpolateBx is not None) and (__InterpolateBy is not None) and
         (__InterpolateBz is not None) ):
        Bx = __InterpolateBx(cartesian_pos)
        By = __InterpolateBy(cartesian_pos)
        Bz = __InterpolateBz(cartesian_pos)
        return np.array([ Bx, By, Bz ]).flatten()
    else:
        meshes = precompute()
        x  = meshes['x']
        y  = meshes['y']
        z  = meshes['z']
        BX = meshes['BX']
        BY = meshes['BY']
        BZ = meshes['BZ']
        
        __InterpolateBx = interpolate.RegularGridInterpolator((x,y,z), BX, bounds_error=False, fill_value=0) # [Tesla]
        __InterpolateBy = interpolate.RegularGridInterpolator((x,y,z), BY, bounds_error=False, fill_value=0) # [Tesla]
        __InterpolateBz = interpolate.RegularGridInterpolator((x,y,z), BZ, bounds_error=False, fill_value=0) # [Tesla]   
        
        return cartesianTesla(cartesian_pos)
