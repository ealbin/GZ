import os
import numpy as np

from gz import coordinates
from gz import units

from gz.results import Results

class HailMary:

    def __init__(self, obj, save_dir='./hailmary'):
        self.filename = save_dir + '/' + os.path.basename(obj.directory.rstrip('/')) + '.txt'
        self.results = obj.results
        atol=1e-2

        p_list  = []
        dp_list = []
        n_list  = []
        dn_list = []
        for r in self.results:
            r.findLastPos()
            p, dp, n, dn = r.getSummary()
            p_list.append(p)
            dp_list.append(dp)
            n_list.append(n)
            dn_list.append(dn)

        p_near  = np.isclose(p_list,  [1.], atol=atol)
        dp_near = np.isclose(dp_list, [1.], atol=atol)
        n_near  = np.isclose(n_list,  [1.], atol=atol)
        dn_near = np.isclose(dn_list, [1.], atol=atol)

        p_dp_both = p_near * dp_near
        n_dn_both = n_near * dn_near

        p_xor_dp = (p_near * ~dp_near) + (~p_near * dp_near)
        n_xor_dn = (n_near * ~dn_near) + (~n_near * dn_near)

        p_neither = ~(p_dp_both + p_xor_dp)
        n_neither = ~(n_dn_both + n_xor_dn)

        p_both = len(p_dp_both[p_dp_both])
        n_both = len(n_dn_both[n_dn_both])
        p_solo = len(p_xor_dp[p_xor_dp])
        n_solo = len(n_xor_dn[n_xor_dn])
        p_none = len(p_neither[p_neither])
        n_none = len(n_neither[n_neither])

        with open(self.filename, 'w') as f:
            f.write('N sims: ' + str(len(self.results)) + "\n")
            f.write('Proton both:    ' + str(p_both) + ', ' + str(p_both / len(self.results) * 100.) + '%' + "\n")
            f.write('One, not both:  ' + str(p_solo) + ', ' + str(p_solo / len(self.results) * 100.) + '%' + "\n")
            f.write('Proton none:    ' + str(p_none) + ', ' + str(p_none / len(self.results) * 100.) + '%' + "\n")
            f.write("\n")
            f.write('Neutron both:   ' + str(n_both) + ', ' + str(n_both / len(self.results) * 100.) + '%' + "\n")
            f.write('One, not both:  ' + str(n_solo) + ', ' + str(n_solo / len(self.results) * 100.) + '%' + "\n")
            f.write('Neutron none:   ' + str(n_none) + ', ' + str(n_none / len(self.results) * 100.) + '%' + "\n")
            f.write('\n')
            
            p_dp_dist = []
            p_record = []
            for pair in np.asarray(self.results)[p_dp_both]:
                _p  = pair.p_last  - coordinates.Cartesian.earth
                _dp = pair.dp_last - coordinates.Cartesian.earth
                
                #!!!! KILOMETERS !!!!
                _dist = Results.HaversineSeparation(_p, _dp)
                _dist = _dist * units.SI.radius_earth / 1000.
                p_dp_dist.append(_dist)
                
                _p  = _p  / np.sqrt(np.dot(_p,  _p ))
                _dp = _dp / np.sqrt(np.dot(_dp, _dp))
                _p_th  = np.rad2deg(np.arccos(_p[2] ))
                _dp_th = np.rad2deg(np.arccos(_dp[2]))
                _p_ph  = np.rad2deg(np.arctan2(_p[1],  _p[0] ))
                _dp_ph = np.rad2deg(np.arctan2(_dp[1], _dp[0]))
                p_record.append([_p_th,  _p_ph,  _dist])
                p_record.append([_dp_th, _dp_ph, _dist])                

            n_dn_dist = []
            n_record = []
            for pair in np.asarray(self.results)[n_dn_both]:
                _n  = pair.n_last  - coordinates.Cartesian.earth
                _dn = pair.dn_last - coordinates.Cartesian.earth
                
                #!!!! KILOMETERS !!!!
                _dist = Results.HaversineSeparation(_n, _dn)
                _dist = _dist * units.SI.radius_earth / 1000.
                n_dn_dist.append(_dist)

                _n  = _n  / np.sqrt(np.dot(_n,  _n ))
                _dn = _dn / np.sqrt(np.dot(_dn, _dn))
                _n_th  = np.rad2deg(np.arccos(_n[2] ))
                _dn_th = np.rad2deg(np.arccos(_dn[2]))
                _n_ph  = np.rad2deg(np.arctan2(_n[1],  _n[0] ))
                _dn_ph = np.rad2deg(np.arctan2(_dn[1], _dn[0]))
                n_record.append([_n_th,  _n_ph,  _dist])
                n_record.append([_dn_th, _dn_ph, _dist])                
                                   
            mean = np.mean(p_dp_dist)
            maxx = np.max(p_dp_dist)
            minn = np.min(p_dp_dist)
            f.write('Proton pair mean, max, min: ' + "\n")
            f.write('\t' + str(mean) + "\n")
            f.write('\t' + str(maxx) + "\n")
            f.write('\t' + str(minn) + "\n")
            f.write("\n")

            mean = np.mean(n_dn_dist)
            maxx = np.max(n_dn_dist)
            minn = np.min(n_dn_dist)
            f.write('Neutron pair mean, max, min: ' + "\n")
            f.write('\t' + str(mean) + "\n")
            f.write('\t' + str(maxx) + "\n")
            f.write('\t' + str(minn) + "\n")
            f.write('\n\n')
            
            f.write('# Proton / Z-1 Fragment\n')
            f.write('# theta [deg], phi [deg], separation [km]\n')
            for th, ph, d in p_record:
                f.write('{} {} {}\n'.format(th, ph, d))
            f.write('\n\n')

            f.write('# Neutron / Z Fragment\n')
            f.write('# theta [deg], phi [deg], separation [km]\n')
            for th, ph, d in n_record:
                f.write('{} {} {}\n'.format(th, ph, d))
            