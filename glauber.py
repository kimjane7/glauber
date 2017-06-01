from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy as np

def get_array(nmax, filename):
	Z = []
	file = open(filename, 'r')
	for line in file.readlines():
		if len(Z)<nmax:
			Z.append([])
			for i in line.split():
				if len(Z[-1])<nmax:
					Z[-1].append(float(i))
	Z = np.array(Z)
	return Z

# set up figure
plt.figure(1)
plt.figure(figsize=(8,4))
plt.rc('text', usetex=True)

# plot contour map of thickness as function of position

T_data = get_array(100,"T_test.dat")
T = plt.subplot(121)
plt.contourf(T_data,100,cmap=plt.cm.get_cmap('magma_r'))
plt.title("Glauber Thickness of Gold Nucleus", color='k', fontsize=20)
plt.contour(T_data,5,colors='#FFFFFF')
T.get_xaxis().set_visible(False)
T.get_yaxis().set_visible(False)
#plt.clabel(T_lines, inline=1, fontsize=10, manual=True)


# plot contour map of energy density of Au-Pb collision pair
b = 0.5
R_Au = 7.056039
R_Pb = 7.169263
x = np.arange(-1.25*R_Au-0.5*b,1.25*R_Pb+0.5*b,(1.25*(R_Au+R_Pb)+b)/300)
y = np.arange(-1.25*R_Pb,1.25*R_Pb,2.5*R_Pb/300)
X,Y = np.meshgrid(x,y)
eps_data = get_array(300,"eps_test.dat")
eps = plt.subplot(122)
plt.contourf(X,Y,eps_data,100,cmap=plt.cm.get_cmap('magma_r'))
plt.title("Energy Density Au-Pb", color='k', fontsize=20)
plt.contour(X,Y,eps_data,5,colors='#FFFFFF')
eps.get_xaxis().set_visible(False)
eps.get_yaxis().set_visible(False)
#plt.clabel(eps_lines, inline=1, fontsize=10, manual=True)

plt.show()