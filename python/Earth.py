#!/usr/bin/env python

"""blablabla
"""

import numpy as np
import os
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime

import Constants
import File

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Development"


def surface_coordinates(TX, PX, R, Zparent, E, proton=False, simulation_dir='../sim/'):
    """ Compute the location (x,y,z) on Earth of where a simulated trajectory ended.
    Earth is taken to be (0,0,0) with the same cartesian directions defined for the solar coordinate
    system (aka, z is "up", x is directed from sun-to-earth-to-space, and y is in the direction of 
    Earth's orbit.
    Simulation directory structure assumed to be:
        simulation_dir/THETAX_???/Z??_E???e18_PHIX???_R?.?.[dat or pdat]
        
    Arguments to identify the simulation file:
        TX      === THETAX, the angle w.r.t. the x-axis [degrees]
        PX      === PHIX, the azimuthal angle in the z-y plane w.r.t. the y axis [degrees]
        R       === distance from Earth when the simulation began [AU]
        Zparent === the parent nucleus atomic number prior to dissintegrating [unit-less]
        E       === the initial energy of the parent nucleus [electronVolts]
        proton  === True/False to identify if this a proton (True) or daughter nucleus (False)
        simulation_dir === specify a root simulation directory
    
    Returns:
        [x, y, z] cartesian location on the surface of the Earth in [meters]
    """    
    
    # check for THETAX directory
    if not os.path.exists( os.path.join(simulation_dir, 'THETAX_{:03}'.format(TX)) ):
        print 'incorrect THETAX value, please select from: ',
        print [x.strip('THETAX_') for x in os.listdir(os.curdir) if x.find('THETAX') == 0]
        sys.exit(1)
        
    # proton or nucleus    
    ext = '.dat'
    if proton:
        ext = '.pdat'
    import pdb
    # find the file or complain
    dirname  = os.path.join( simulation_dir, 'THETAX_{:03}'.format(TX) )
    filename = 'Z{:02}_E{:03.0f}e18_PHIX{:03.0f}_R{:2.1f}'.format(Zparent, E, PX, R) + ext
    filepath = os.path.join(dirname, filename)
    if not os.path.exists(filepath):
        data = {'Zparent':[], 'E':[], 'PX':[], 'R':[]}
        contents = [ [x[0], x[1], x[2]] for x in os.walk(dirname) ]
        for file in contents[0][2]:
            tokens = file.split('_')
            data[ 'Zparent'].append(   int( tokens[0].strip('Z') ) )
            data[       'E'].append(   int( tokens[1].strip('E').strip('e18') ) )
            data[      'PX'].append(   int( tokens[2].strip('PHIX') ) )
            data[       'R'].append( float( tokens[3].strip('R.pdat') ) )
        print 'incorrect value(s), please select from: TX, Zparent, E, PX, R, proton=False(default)'
        print '\t TX = {}'.format(TX)
        for key in data.keys():
            tag = ''
            if key == 'E':
                tag = 'e18'
            shortlist = sorted( list(set( data[key] )) )
            print '\t ' + key + ' = {}'.format(shortlist) + tag
        sys.exit(1)

    # if all is good, go ahead and fetch from file
    #---------------------------------------------
    earth_pos = [1,0,0] # solar cartesian system [astronomical units]
    file_data = File.read(filepath)
    if file_data['exit_info'] != 'earth-plane':
        print 'simulation did not end at the earth-plane'
        sys.exit(1)
    
    R_earth = Constants.RadiusMeters('earth')        
    x = Constants.AU2meters( file_data['last_pos'] - earth_pos ) # earth cartesian system [meters]
    proximity_ratio = np.sqrt(np.dot(x,x)) / R_earth
    if proximity_ratio > 1.0:
        print 'trajectory misses earth by {:3.1f} x Re'.format(proximity_ratio)
        sys.exit(1)
        
    beta = file_data['last_beta'] 
    beta = beta / np.sqrt( np.dot(beta, beta) )
    x_dot_beta = np.dot(x, beta)
    # solve for scale factor 's' to bring x from the plane to the surface:
    s = x_dot_beta + np.sqrt( R_earth**2. + x_dot_beta**2. - np.dot(x,x) ) 
    return (x - s * beta) # (x,y,z) location on the surface of the Earth [meters]

def latlon(position, lat_offset = 0., lon_offset = 0.):
    """ Returns the latitude and longitude [lat, lon] in degrees for cartesian
    coordinates 'position'.
    """
    pos   = np.array(position)
    p_hat = pos / np.sqrt( np.dot(pos, pos) )
    
    lat = ( 180. / np.pi ) * np.arcsin(  p_hat[2] )
    lon = ( 180. / np.pi ) * np.arctan2( p_hat[1], p_hat[0] )
    
    return np.array([lat + lat_offset, lon + lon_offset])

def plot(latlon):
    """ Draws the Earth and plots a point latlon = [latitude, longitude] in degrees
    """
    fig = plt.figure(figsize=(8, 8))

    lat, lon = latlon

    map = Basemap(projection='ortho',lat_0=39.8283, lon_0=-98.5795)
    map.drawparallels( np.arange(-90, 90, 30), color='w' )
    map.drawmeridians( np.arange(map.lonmin, map.lonmax+30, 60), color='w' )
    map.bluemarble()
    map.drawcoastlines(color='w')

    # shade the night areas, with alpha transparency so the
    # map shows through. Use current time in UTC.
    #date = datetime.utcnow()
    #birthday
    date = datetime(1984, 9, 20, 23, 19, 0, 0, pytz.timezone('US/Pacific'))
    date = date.astimezone(pytz.timezone('UTC'))
    CS=map.nightshade(date, alpha=0.4)
    plt.title('Day/Night Map for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
    
    plt.show()
    