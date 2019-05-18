import propagate
import time
import numpy as np
import matplotlib.pyplot as plt

o = propagate.Outgoing([1,0,0],[-1,0,0],92,238,2e18)

start = time.time()
for _ in range(10000):
    o.propogate()
print(time.time() - start)

t = np.asarray(o.telementry)

plt.plot(t.T[0], t.T[1])
