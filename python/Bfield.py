"""Precompute the total magnetic field, store to disk 
and use it as an interpolated look-up table to profoundly accelerate
numeric integration.
"""
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import os
import SolarMagneticModel
import Transform

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


def cartesianTesla(spacelimit=6, resolution=1000, autoload=True, directory='tables', b_fname='cartesianBfield.Tesla'):
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
    b_path    = os.path.abspath( os.path.join( base_dir, directory, b_fname ) )
    b_exists  = os.path.isfile(b_path)
    if (not b_exists):
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
        for i, (ix, iy, iz) in enumerate( zip( X.flatten(), Y.flatten(), Z.flatten() ) ):
            b_solar = SolarMagneticModel.sumBfieldTesla( np.array([ ix, iy, iz ]) )
            bx[i] = b_solar[0] # [Tesla]
            by[i] = b_solar[1] # [Tesla]
            bz[i] = b_solar[2] # [Tesla]

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


def interpolate( cartesian_pos ):
    """Returns cartesian [Tesla] values (Bx, By, Bz) for cartesian_pos = (x, y, z)
    """
    meshes = cartesianTesla()
    x  = meshes['x']
    y  = meshes['y']
    z  = meshes['z']
    BX = meshes['BX']
    BY = meshes['BY']
    BZ = meshes['BZ']
    
    Bx = RegularGridInterpolator((x,y,z), BX) # [Tesla]
    By = RegularGridInterpolator((x,y,z), BY) # [Tesla]
    Bz = RegularGridInterpolator((x,y,z), BZ) # [Tesla]    

    return np.array([ Bx(cartesian_pos), By(cartesian_pos), Bz(cartesian_pos) ])
