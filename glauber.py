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

plt.contourf(Z,100,cmap=plt.cm.get_cmap('YlGnBu_r'))
plt.contour(Z,5,colors='#FFFFFF')