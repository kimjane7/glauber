from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy as np

def get_array(filename):
	Z = []
	file = open(filename, 'r')
	for line in file.readlines():
		if len(Z)<250:
			Z.append([])
			for i in line.split():
				if len(Z[-1])<250:
					Z[-1].append(float(i))
	Z = np.array(Z)
	return Z

# set up figure
plt.figure(1)
plt.figure(figsize=(8,8))
plt.rc('text', usetex=True)

# plot contour map
Z = get_array("thickness_test.dat")
ax = plt.subplot(111)
plt.contourf(Z,100,cmap=plt.cm.get_cmap('magma_r'))
fig = plt.contour(Z,5,colors='#FFFFFF')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.clabel(fig, inline=1, fontsize=10, manual = True)


plt.show()