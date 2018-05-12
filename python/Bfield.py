"""Precompute the total magnetic field, store to disk 
and use it as an interpolated look-up table to profoundly accelerate
numeric integration.
"""
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import os
import SolarMagneticModel
import Transform

__x  = None
__y  = None
__z  = None
__BX = None
__BY = None
__BZ = None

def cartesianTesla(spacelimit=10, resolution=10, autoload=True, 
                   directory='tables', b_fname='cartesianBfield.Tesla'):
    """Returns total magnetic field X, Y, Z, BX, BY, BZ meshes by
    by disk-read or re-generation.  Field density in Teslas.
    
    spacelimit : radial reach (r) of the space volume 
                (x, y, z) === (-r to r) by (-r to r) by (-r to r) [astronomical units]

    resolution : the number of samples taken between (-r to r) along each dimension.
    
    autoload   : if True, look FIRST to disk for existing spacelimit/resolution.
                    if (no preexisting) or (incorrect spacelimit or resolution):
                        regenerate from scratch and overwrite existing.
                 if False, regenerate from scratch and overwrite existing.
    
    directory  : subdirectory with magnetic field text file
    
    b_fname    : filename for magnetic field text file
    
    returns dictionary { 'x', 'y', 'z', 'BX', 'BY', 'BZ' }

    """
    spacelimit = int( spacelimit ) # [astronomical units] (integer for easy file read)
    resolution = int( resolution ) # N divisions (integer for easy file read)

    x = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    y = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]
    z = np.linspace(-spacelimit, spacelimit, resolution) # [astronomical units]

    spacelimit = np.array([spacelimit])
    resolution = np.array([resolution])

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    shape = np.array(X.shape)
    size = X.flatten().size
    
    bx = np.zeros(size) # [Tesla]
    by = np.zeros(size) # [Tesla]
    bz = np.zeros(size) # [Tesla]    

    regen = False
    if not autoload:
        regen = True

    directory = 'tables'
    b_fname   = 'cartesianBfield.Tesla'
    b_path    = os.path.abspath( os.path.join( os.path.curdir, directory, b_fname ) )
    b_exists  = os.path.isfile(b_path)
    if (not b_exists):
        regen = True

    text_sep = ', '
    
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

    # for whatever circumstance, regenerate and overwrite
    if regen:
        for i, (ix, iy, iz) in enumerate( zip( X.flatten(), Y.flatten(), Z.flatten() ) ):
            b_solar = SolarMagneticModel.sumBfieldTesla( np.array([ ix, iy, iz ]) )

            bx[i] = b_solar[0] # [Tesla]
            by[i] = b_solar[1] # [Tesla]
            bz[i] = b_solar[2] # [Tesla]

        with open(b_path, 'w') as b_f:
            header = ( 'rows: 0:this header, 1:spacelimit [AU], 2:resolution, 3:shape, '
                       '4:x [AU], 5:y [AU], 6:z [AU], 7:BX [T], 8:BY [T], 9:BZ [T]' )
            b_f.write(header + '\n')

            spacelimit.tofile(b_f, sep=text_sep)
            b_f.write('\n')

            resolution.tofile(b_f, sep=text_sep)
            b_f.write('\n')

            shape.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            x.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            y.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            z.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            bx.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            by.tofile(b_f, sep=text_sep)
            b_f.write('\n')
            
            bz.tofile(b_f, sep=text_sep)
            b_f.write('\n')
        
    BX = bx.reshape(shape)
    BY = by.reshape(shape)
    BZ = bz.reshape(shape)
   
    global __x, __y, __z, __BX, __BY, __BZ
    __x  = x
    __y  = y
    __z  = z
    __BX = BX
    __BY = BY
    __BZ = BZ
    
    return {'x':x, 'y':y, 'z':z, 'BX':BX, 'BY':BY, 'BZ':BZ}    


def interpolate( cartesian_pos ):
    """Returns cartesian [Tesla] values (Bx, By, Bz) for cartesian_pos = (x, y, z)
    """
    global __x, __y, __z, __BX, __BY, __BZ       
    if ( (__x is None )  or (__y is None)  or (__z is None) or
         (__BX is None ) or (__BY is None) or (__BZ is None) ):
       meshes = cartesianTesla()
       __x  = meshes['x']
       __y  = meshes['y']
       __z  = meshes['z']
       __BX = meshes['BX']
       __BY = meshes['BY']
       __BZ = meshes['BZ']
    
    Bx = RegularGridInterpolator(cartesian_pos, __BX) # [Tesla]
    By = RegularGridInterpolator(cartesian_pos, __BY) # [Tesla]
    Bz = RegularGridInterpolator(cartesian_pos, __BZ) # [Tesla]    

    return np.array([ Bx, By, Bz ])
