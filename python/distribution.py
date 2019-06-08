
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

class D:

    def __init__(self, filename):

        with open(filename) as file:
            
            self.pos = []
            self.beta = []
            for line in file.readlines():
                tokens = line.split()
                self.pos.append(np.asarray(tokens[:3], dtype=np.float64))
                self.beta.append(np.asarray(tokens[3:6], dtype=np.float64))
            
    def getSolar(self):
        self.pos_r = []
        R = 1000.
        for pb in zip(self.pos, self.beta):
            pdotb = np.dot(pb[0], pb[1])
            p = np.sqrt(np.dot(pb[0], pb[0]))
            r = -pdotb + np.sqrt(R*R - p*p + pdotb*pdotb)
            
            a = pb[0] + r * pb[1]
            d = np.sqrt(np.dot(a,a))
            if (not np.isclose([d], [R])):
                print("Not close: " + str(d))
            else:
                self.pos_r.append(a)
        
        
    def getGalactic(self):
        self.pos_g = []
        R = 1e12
        mw = np.asarray([1.65e9, 0, 0])
        for pb in zip(self.pos, self.beta):
            pp = mw + pb[0]
            pdotb = np.dot(pp, pb[1])
            p = np.sqrt(np.dot(pp, pp))
            r = -pdotb + np.sqrt(R*R - p*p + pdotb*pdotb)
            
            a = pp + r * pb[1]
            d = np.sqrt(np.dot(a,a))
            if (not np.isclose([d], [R])):
                print("Not close: " + str(d))
            else:
                self.pos_g.append(a)
        
            
    def plotSolar(self, weighting=False):
        lat = []
        lon = []
        for pos in self.pos_r:
            p = np.sqrt(np.dot(pos, pos))
            hat = pos / p        
            theta = 180. / np.pi * np.arccos(hat[2])
            phi = 180. / np.pi * np.arctan2(hat[1], hat[0])
            
            lat.append(90. - theta)
            if (phi <= 180.):
                lon.append(phi)
            else:
                lon.append(phi - 360.)
        
        weights = []
        if (weighting):
                for _ in zip(lon, lat):
                    ln = _[0]
                    lt = _[1]
                    theta = 90 - lt
                    weights.append(1. / np.sin(np.deg2rad(theta)))
        else:
            weights = None
        
        fig1 = plt.figure(figsize=[20,20])
        lon_bins = np.concatenate([np.linspace(-180,0,181), np.linspace(1,180,180)])
        lat_bins = np.concatenate([np.linspace(-90, 90, 181)])
        plt.hist2d(lon, lat, bins=[lon_bins, lat_bins], normed=False, cmin=0., weights=weights, vmin=0, vmax=500)
        plt.colorbar()
        plt.title('Solar Sky: 1000 AU from the Sun')
        plt.xlabel('Solar Longitude [deg]')
        plt.ylabel('Solar Latitude [deg]')
        
    def plotGalactic(self):
        lat = []
        lon = []
        for pos in self.pos_g:
            p = np.sqrt(np.dot(pos, pos))
            hat = pos / p        
            theta = 180. / np.pi * np.arccos(hat[2])
            phi = 180. / np.pi * np.arctan2(hat[1], hat[0])
            
            lat.append(90. - theta)
            if (phi <= 180.):
                lon.append(phi)
            else:
                lon.append(phi - 360.)
        
        fig1 = plt.figure(figsize=[20,20])
        lon_bins = np.concatenate([np.linspace(-180,0,181), np.linspace(1,180,180)])
        lat_bins = np.concatenate([np.linspace(-90, 90, 181)])
        plt.hist2d(lon, lat, bins=[lon_bins, lat_bins], normed=False, cmin=0.)
        plt.colorbar()
        plt.title('Galactic Sky: 1e12 AU from the Sun')
        plt.xlabel('Galactic Longitude [deg]')
        plt.ylabel('Galactic Latitude [deg]')
        