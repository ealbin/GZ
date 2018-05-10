from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import SolarMagneticModel

## mesh
x = np.linspace(-5, 5, 20)
y = np.linspace(-5, 5, 20)
z = np.linspace(-5, 5, 20)

X, Y, Z = np.meshgrid(x, y, z)
x = X.reshape(X.size)
y = Y.reshape(Y.size)
z = Z.reshape(Z.size)

Bx = []
By = []
Bz = []
for i in range(x.size):
    B = SolarMagneticModel.solarBfield(np.array([ x[i], y[i], 0 ]))
    Bx.append(B[0])
    By.append(B[1])
    Bz.append(B[2])
    
Bx = np.array(Bx).reshape(X.shape)
By = np.array(By).reshape(Y.shape)
Bz = np.array(Bz).reshape(Z.shape)    
    
## figure
fig = plt.figure()
ax  = fig.gca(projection='3d')
ax.quiver(X, Y, Z, Bx, By, Bz, color='b', length=1, normalize=True)
plt.show()
