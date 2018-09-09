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
import PDFfield

__author__ = "Eric Albin"
__copyright__ = "Copyright 2018, The CRAYFIS Project"
__credits__ = ["Eric Albin"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Eric Albin"
__email__ = "Eric.K.Albin@gmail.com"
__status__ = "Development"


def surface_coordinates(filepath):
    """ Compute the location (x,y,z) on Earth of where a simulated trajectory ended.
    Earth is taken to be (0,0,0) with the same cartesian directions defined for the solar coordinate
    system (aka, z is "up", x is directed from sun-to-earth-to-space, and y is in the direction of 
    Earth's orbit.
    Parameters:
        filepath === a string filepath for a simulation result
    Returns:
        [x, y, z] cartesian location on the surface of the Earth in [meters]
    """    
    # check path
    if not os.path.exists( filepath ):
        print 'filepath {} does not exist'.format(filepath)
        return None
        
    earth_pos = [1,0,0] # solar cartesian system [astronomical units]
    file_data = File.read(filepath)
    if file_data['exit_info'] != 'earth-plane':
        print 'simulation {} did not end at the earth-plane'.format(filepath)
        return None
    
    R_earth = Constants.RadiusMeters('earth')        
    x = Constants.AU2meters( file_data['last_pos'] - earth_pos ) # earth cartesian system [meters]
    proximity_ratio = np.sqrt(np.dot(x,x)) / R_earth
    if proximity_ratio > 1.0:
        print 'trajectory {} misses earth by {:3.1f} x Re'.format(filepath, proximity_ratio)
        return None
        
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
    
    
def fetch_sims(simulation_dir='../sim/'):
    """ Returns a 'master list' of all simulation files in the simulation_dir.
    """
    master_list = []
    for root, dirs, files in os.walk(simulation_dir):
        if root.find('THETAX_') < 0:
            continue
        if len(dirs) > 0:
            continue

        for file in files:
            # catch and remove protons from the list
            if file.find('.pdat') > -1:
                continue
            
            # fill in master_list format
            data = {'TX':None, 'PX':None, 'R':None, 'Zparent':None, 'Aparent':None, 'Eparent':None, 
                    'Edaughter':None, 'Enucleon':None, 'path':None,
                    'lat':None, 'lon':None, 'lat0':None, 'lon0':None, 'lat1':None, 'lon1':None}    
            data['TX']   = int(root[-3:])
            data['path'] = os.path.abspath( os.path.join( root, file ) )
            
            tokens = file.split('_')
            data['Zparent'] = int(   tokens[0].strip('Z') )
            data['Eparent'] = int(   tokens[1].strip('E').strip('e18') )
            data[     'PX'] = int(   tokens[2].strip('PHIX') )
            data[      'R'] = float( tokens[3].strip('R.pdat') )
            
            A = Constants.mass_number(data['Zparent'])
            data['Aparent'] = A
            data['Edaughter'] = (A - 1.) / float(A) * data['Eparent']
            data['Enucleon']  =  1. / float(A) * data['Eparent']
            
            position = surface_coordinates(data['path'])
            if position is None:
                print '\tskipping..'
                continue
            lat, lon = latlon(position) 
            data['lat'] = lat
            data['lon'] = lon

            # neutron ejection information
            R_earth = Constants.RadiusMeters('earth')        
            tx = data['TX'] * np.pi / 180.
            px = data['PX'] * np.pi / 180.
            x  = R_earth * np.cos(tx)
            y  = R_earth * np.sin(tx) * np.cos(px)
            z  = R_earth * np.sin(tx) * np.sin(px)
            lat0, lon0 = latlon([x, y, z])
            data['lat0'] = lat0
            data['lon0'] = lon0
        
            # proton ejection information
            position = surface_coordinates(data['path'][:-3] + 'pdat')
            if position is None:
                print '\t skipping b/c proton..'
                continue
            lat, lon = latlon(position) 
            data['lat1'] = lat
            data['lon1'] = lon
        
            master_list.append(data)

    return master_list


def filter_sims(master_list, TX=None, PX=None, R=None, Zparent=None, Eparent=None):
    """ Filters master_list for parameters specified, e.g.:
    TX=5, R=[1,2,3] ==> return only simulations with THETAX=5 and R=(1 or 2 or 3)
    """
    def ismatch(search_this, match_this):
        values = np.array(match_this).flatten()
        if search_this in values:
            return True
        return False
        
    keep_list = []
    for item in master_list:
        if (TX is not None) and (not ismatch(item['TX'], TX)):
            continue
        if (PX is not None) and (not ismatch(item['PX'], PX)):
            continue
        if ( R is not None) and (not ismatch(item[ 'R'],  R)):
            continue
        if (Eparent is not None) and (not ismatch(item['Eparent'], Eparent)):
            continue
        if (Zparent is not None) and (not ismatch(item['Zparent'], Zparent)):
            continue
        keep_list.append(item)           
    
    return keep_list


def add_probabilities(listing):
    for item in listing:
        TX = item['TX']             # [degrees]
        PX = item['PX']             # [degrees]
        R  = item['R']              # [AU]
        A  = item['Aparent']        # [unitless]
        E  = item['Eparent'] * 1e18 # [electronVolts]
        item['probability'] = PDFfield.integratePath(TX, PX, R, A, E)
    
    
def HaversineNormalized(lat1, lon1, lat2, lon2):
    """ Normalized the earth radius aka 1 = Re
    """
    lat1 *= np.pi / 180.
    lon1 *= np.pi / 180.
    lat2 *= np.pi / 180.
    lon2 *= np.pi / 180.
    
    part1 = np.sin( (lat2 - lat1) / 2. )**2.
    part2 = np.cos(lat1) * np.cos(lat2)
    part3 = np.sin( (lon2 - lon1) / 2. )**2.

    return 2. * np.arcsin( np.sqrt(part1 + part2 * part3) )
    

def plot_histogram(listing):
    """
    """
    fig  = plt.figure(figsize=(8, 8))
    bins = np.logspace(-3, 0, 50) 
    dist_p = []
    dist_n = []
    weight = []
    for entity in listing:
        lat  = entity['lat']
        lon  = entity['lon']
        lat0 = entity['lat0']
        lon0 = entity['lon0']
        lat1 = entity['lat1']
        lon1 = entity['lon1']
        prob = entity['probability']['probability']
        
        weight.append(prob)
        
        # distance between daughter and proton
        dist_p.append( HaversineNormalized( lat, lon, lat1, lon1 ) )
        
        # distance between daughter and neutron
        dist_n.append( HaversineNormalized( lat, lon, lat0, lon0 ) )
    
    plt.hist( dist_p, bins=bins, normed=True, weights=weight, fc=(0, 0, 1, .5) )
    plt.hist( dist_n, bins=bins, normed=True, weights=weight, fc=(1, 0, 0, .5) )
    plt.gca().set_xscale('log')
    plt.axvline(x=np.sqrt(3000)*1000. / Constants.RadiusMeters('earth'))
    plt.show()
    
def plot_triangle(plot_list):
    """
    """
    fig = plt.figure(figsize=(10, 10))
    for entity in plot_list:
        lat  = entity['lat']
        lon  = entity['lon']
        lat0 = entity['lat0']
        lon0 = entity['lon0']
        lat1 = entity['lat1']
        lon1 = entity['lon1']
        
        # distance between daughter and proton
        dist_p = HaversineNormalized( lat, lon, lat1, lon1 )
        
        # distance between daughter and neutron
        dist_n = HaversineNormalized( lat, lon, lat0, lon0 )
        
        angle = (entity['R'] / 6.) * 90. * np.pi / 180.
        x_p = dist_p * np.cos(angle)
        y_p = dist_p * np.sin(angle)
        x_n = dist_n * np.sin(angle-np.pi/2.)
        y_n = dist_n * np.cos(angle-np.pi/2.)
        
        plt.xlim([-1,1])
        plt.ylim([0,1])
        
        #x = [ 0, dist_p, 0, 0 ]
        #y = [ 0, 0, dist_n, 0 ]
        #plt.plot(x, y)
        
        plt.scatter([0], [0], s=2, marker='o', color='g')
        plt.scatter(x_p, y_p, marker='X', color='b')
        plt.scatter(x_n, y_n, marker='D', color='r')
        
    plt.show()

def plot(plot_list):
    """ Draws the Earth and plots a point latlon = [latitude, longitude] in degrees
    """
    fig = plt.figure(figsize=(20, 10))

    map = Basemap(projection='hammer', lon_0=180)
    map.bluemarble()    
    map.drawparallels( np.arange(-90, 90, 30), color='w' )
    map.drawmeridians( np.arange(map.lonmin, map.lonmax+30, 60), color='w' )
    
    # plot key:
    neutron = 'm'
    proton  = 'c'
    helium  = 'y'
    oxygen  = 'b'
    iron    = 'r'
    uranium = 'g'
    
    r1 = '^' # triangle up
    r2 = 's' # square
    r3 = 'p' # pentagon
    r4 = 'h' # hexagon 1
    r5 = '8' # octagon
    r6 = 'o' # circle
    
    for item in plot_list:
        color  = None
        if item['proton'] == True:
            color = proton
        elif item['Zparent'] == 4:
            color = helium
        elif item['Zparent'] == 8:
            color = oxygen
        elif item['Zparent'] == 26:
            color = iron
        elif item['Zparent'] == 92:
            color = uranium
        else:
            print 'uh oh, weird Z'
            return

        marker = None        
        if item['R'] <= 1.0:
            marker = r1
        elif item['R'] <= 2.0:
            marker = r2
        elif item['R'] <= 3.0:
            marker = r3
        elif item['R'] <= 4.0:
            marker = r4
        elif item['R'] <= 5.0:
            marker = r5
        elif item['R'] > 5.0:
            marker = r6
        else:
            print 'uh oh, weird R'
            return
        
        x, y = map(item['lon'], item['lat'])
        map.plot(x, y, marker=marker, color=color)
        
        x0, y0 = map(item['lon0'], item['lat0'])
        map.plot(x0, y0, marker='x', color=neutron)
        
    plt.show()    

    #map.drawmapboundary(fill_color='aqua')
    #map.fillcontinents(color='coral',lake_color='aqua')
    #map.drawcoastlines()
    #x, y = map(lon, lat)
    #map.plot(x,y, marker='D',color='m')
    
    #map = Basemap(projection='ortho',lat_0=39.8283, lon_0=-98.5795)
    #map.drawcoastlines(color='w')

    # shade the night areas, with alpha transparency so the
    # map shows through. Use current time in UTC.
    #date = datetime.utcnow()
    #birthday
    #date = datetime(1984, 9, 20, 23, 19, 0, 0, pytz.timezone('US/Pacific'))
    #date = date.astimezone(pytz.timezone('UTC'))
    #CS=map.nightshade(date, alpha=0.4)
    #plt.title('Day/Night Map for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
    
















    #def surface_coordinates(TX, PX, R, Zparent, E, proton=False, simulation_dir='../sim/'):
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
    """
