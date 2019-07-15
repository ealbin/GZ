

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import gz

def why_sine(N=10000):
    
    def xyz(theta, phi):
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
        return x, y, z
    
    tu = np.pi * np.random.random(size=N)
    pu = 2*np.pi * np.random.random(size=N)

    ts, ps = gz.earth.Earth.randomThetaPhi(N)

    plt.rc('font', size=20)
    fig = plt.figure(figsize=[24,12])

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    x, y, z = xyz(tu, pu)
    ax.scatter3D(x, y, z)
    plt.title('Uniform')
    ax.view_init(45, 0)
    
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    x, y, z = xyz(ts, ps)
    ax.scatter3D(x, y, z)
    plt.title('Sin(theta) weighting')
    ax.view_init(45, 0)
    
    plt.tight_layout()
    
    
def trajectory(directory, scale=10):
    
    def rel_earth(telemetry):
        earth = np.asarray([1., 0., 0., 0., 0., 0., 0.])
        x, y, z, bx, by, bz, _ = (telemetry - earth).T
        return x, y, z
    
    def local_xyz(tel1, tel2):
        # treates tel1 as x
        
        x1, y1, z1, bx1, by1, bz1, _ = tel1[-1]
        x2, y2, z2, bx2, by2, bz2, _ = tel2[-1]
        
        b1 = np.asarray([bx1, by1, bz1])
        b2 = np.asarray([bx2, by2, bz2])
        
        x = b1 / np.sqrt(np.dot(b1, b1))
        
        z = np.cross(x, b2)
        z = z / np.sqrt(np.dot(z, z))
        
        y = np.cross(z, x)
        y = y / np.sqrt(np.dot(y, y))
        
        return x, y, z

    def get_proj(tx, ty, tz, ux, uy, uz):
        t = np.asarray([tx, ty, tz])
        x = np.dot(ux, t)
        y = np.dot(uy, t)
        z = np.dot(uz, t)
        return x, y, z
    
    R = gz.results.Results(directory=directory, full_telemetry=True)
    er = gz.units.SI.radius_earth * gz.units.Change.meter_to_AU

    total = len(R.results)        
    for count, r in enumerate(R.results):
        
        print('\rplotting.. {:.2f}%       '.format(100.*(count+1)/float(total)), end='\r', flush=True)
        p = r.p_telemetry
        n = r.n_telemetry
        dp = r.dp_telemetry
        dn = r.dn_telemetry
         
        px, py, pz = rel_earth(p)
        nx, ny, nz = rel_earth(n)
        dpx, dpy, dpz = rel_earth(dp)
        dnx, dny, dnz = rel_earth(dn)
         
        ux, uy, uz = local_xyz(p, dp)
        px, py, pz = get_proj(px, py, pz, ux, uy, uz)
        nx, ny, nz = get_proj(nx, ny, nz, ux, uy, uz)
        dpx, dpy, dpz = get_proj(dpx, dpy, dpz, ux, uy, uz)
        dnx, dny, dnz = get_proj(dnx, dny, dnz, ux, uy, uz)
        
        px, py, pz = px/er, py/er, pz/er
        nx, ny, nz = nx/er, ny/er, nz/er
        dpx, dpy, dpz = dpx/er, dpy/er, dpz/er
        dnx, dny, dnz = dnx/er, dny/er, dnz/er
        
        plt.rc('font', size=20)
        fig = plt.figure(figsize=[30,20])
    
        ################################## Top ################################
    
        ax = fig.add_subplot(2, 3, 1)
        order = [ [0., 'earth'],
                  [(pz[0]  + pz[-1] ) / 2., [px,  py ], 'r:', 'proton' ],
                  [(nz[0]  + nz[-1] ) / 2., [nx,  ny ], 'g:', 'neutron' ],
                  [(dpz[0] + dpz[-1]) / 2., [dpx, dpy], 'r', 'Z-1 p-daughter'],
                  [(dnz[0] + dnz[-1]) / 2., [dnx, dny], 'g', 'Z n-daughter'] ]
        order.sort()
        plots = {}
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Top-Down View', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory x [earth radii]')
        ax.set_ylabel('trajectory y [earth radii]')
        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, scale)
        ax.set_aspect('equal')
        ax.legend( (plots['r:'], plots['r'], plots['g:'], plots['g']),
                   ('Proton', 'Z-1 Fragment', 'Neutron', 'Z Fragment') )


        ################################ Bottom ###############################
    
        ax = fig.add_subplot(2, 3, 4)
        order = [ [0., 'earth'],
                  [(pz[0]  + pz[-1] ) / 2., [px,  py ], 'r:', 'proton' ],
                  [(nz[0]  + nz[-1] ) / 2., [nx,  ny ], 'g:', 'neutron' ],
                  [(dpz[0] + dpz[-1]) / 2., [dpx, dpy], 'r', 'Z-1 p-daughter'],
                  [(dnz[0] + dnz[-1]) / 2., [dnx, dny], 'g', 'Z n-daughter'] ]
        order.sort(reverse=True)
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Bottom-Up View', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory x')
        ax.set_ylabel('trajectory y')
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlim(-scale, scale)
        ax.set_ylim(scale, -scale)
        ax.set_aspect('equal')
        
        ################################ Left #################################
    
        ax = fig.add_subplot(2, 3, 2)
        order = [ [0., 'earth'],
                  [(py[0]  + py[-1] ) / 2., [px,  pz ], 'r:', 'proton' ],
                  [(ny[0]  + ny[-1] ) / 2., [nx,  nz ], 'g:', 'neutron' ],
                  [(dpy[0] + dpy[-1]) / 2., [dpx, dpz], 'r', 'Z-1 p-daughter'],
                  [(dny[0] + dny[-1]) / 2., [dnx, dnz], 'g', 'Z n-daughter'] ]
        order.sort(reverse=True)
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Look Left', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory x')
        ax.set_ylabel('trajectory z')
        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, scale)
        ax.set_aspect('equal')

        ################################ Right ################################
    
        ax = fig.add_subplot(2, 3, 5)
        order = [ [0., 'earth'],
                  [(py[0]  + py[-1] ) / 2., [px,  pz ], 'r:', 'proton' ],
                  [(ny[0]  + ny[-1] ) / 2., [nx,  nz ], 'g:', 'neutron' ],
                  [(dpy[0] + dpy[-1]) / 2., [dpx, dpz], 'r', 'Z-1 p-daughter'],
                  [(dny[0] + dny[-1]) / 2., [dnx, dnz], 'g', 'Z n-daughter'] ]
        order.sort(reverse=False)
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Look Right', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory x')
        ax.set_ylabel('trajectory z')
        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, scale)
        ax.set_aspect('equal')

        ################################ Head-On ##############################
    
        ax = fig.add_subplot(2, 3, 3)
        order = [ [0., 'earth'],
                  [(px[0]  + px[-1] ) / 2., [py,  pz ], 'r:', 'proton' ],
                  [(nx[0]  + nx[-1] ) / 2., [ny,  nz ], 'g:', 'neutron' ],
                  [(dpx[0] + dpx[-1]) / 2., [dpy, dpz], 'r', 'Z-1 p-daughter'],
                  [(dnx[0] + dnx[-1]) / 2., [dny, dnz], 'g', 'Z n-daughter'] ]
        order.sort(reverse=False)
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Head-On', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory y')
        ax.set_ylabel('trajectory z')
        ax.set_xlim(-scale, scale)
        ax.set_ylim(-scale, scale)
        ax.set_aspect('equal')
        
        ################################ Tail-Forward #########################
    
        ax = fig.add_subplot(2, 3, 6)
        order = [ [0., 'earth'],
                  [(px[0]  + px[-1] ) / 2., [py,  pz ], 'r:', 'proton' ],
                  [(nx[0]  + nx[-1] ) / 2., [ny,  nz ], 'g:', 'neutron' ],
                  [(dpx[0] + dpx[-1]) / 2., [dpy, dpz], 'r', 'Z-1 p-daughter'],
                  [(dnx[0] + dnx[-1]) / 2., [dny, dnz], 'g', 'Z n-daughter'] ]
        order.sort(reverse=True)
        for zorder, _ in enumerate(order):
            if (_[1] == 'earth'):
                earth = plt.Circle((0,0), 1., color='b', zorder=zorder)
                ax.add_artist(earth)
            else:
                x, y = _[1]
                plots[_[2]], = ax.plot(x, y, _[2], label=_[3], zorder=zorder, linewidth=3.0)

        ax.text(.02, .95, 'Tail-Forward', transform=ax.transAxes, fontweight='bold')
        ax.set_xlabel('trajectory y')
        ax.set_ylabel('trajectory z')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.set_xlim(scale, -scale)
        ax.set_ylim(-scale, scale)
        ax.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(r.dirname + '/plots/' + str(scale) + '_' + r.filename.rstrip('incoming')+'png')
        plt.close()