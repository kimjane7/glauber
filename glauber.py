#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt 
import numpy as np
import os

def get_array(nmax, filename):
	Z = []
	file = open(filename, 'r')
	for line in file.readlines():
		if len(Z) < nmax:
			Z.append([])
			for i in line.split():
				if len(Z[-1])<nmax:
					Z[-1].append(float(i))
	Z = np.array(Z)
	return Z

# set up figure
plt.figure(1)
plt.figure(figsize=(7,7))
plt.rc('text', usetex=True)

# run the glauber program to generate the data file
os.system("./Glauber")

# open file generated by the glauber program
data = get_array(100, "energy_density.dat")

# plot contour map of thickness as function of position
T = plt.subplot(111)
plt.contourf(data, 100, cmap=plt.cm.get_cmap('magma_r'))
plt.contour(data, 5, colors='#FFFFFF')
T.get_xaxis().set_visible(False)
T.get_yaxis().set_visible(False)
#plt.clabel(T_lines, inline=1, fontsize=10, manual=True)

plt.savefig('energy_density.pdf',format='pdf')
os.system('open energy_density.pdf')