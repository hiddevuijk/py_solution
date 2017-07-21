import numpy as np
import matplotlib.pyplot as plt
from sys import exit


gr = np.loadtxt("gr.dat")
g0 = np.loadtxt("g0.dat")
r = np.loadtxt("r.dat")
plt.plot(r,gr,label="g(r)")
#plt.plot(r,g0,label="g0")

plt.xlim([0.7,2.5])

plt.legend()
plt.tight_layout()
plt.show()



