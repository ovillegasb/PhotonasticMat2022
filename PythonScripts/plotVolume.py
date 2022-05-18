
import numpy as np
import matplotlib.pyplot as plt


box = np.loadtxt("box.dat")

# mass from DensPiston script
massTOT = 9.135820996651693e-23
# vol in m3
volume = box[:, 0] * box[:, 1] * box[:, 2] * 1e-30

density = massTOT / volume

plt.plot(density)
plt.show()

