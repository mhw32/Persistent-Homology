# -- Figure 2 : toy model demo -- 

import numpy as np
from scipy.ndimage.filters import gaussian_filter
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.numpy2ri import numpy2ri
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set_style('whitegrid')


def rds_to_np(Rfile):
    ''' Convert .RData to be able to load 
        the array into Python.

        Rfile := (str) location of file to 
        translate to Python.
    ''' 
    raw = robjects.r['readRDS'](Rfile)
    return raw


raw = rds_to_np('intermediate/fig_2_tmp.rds')
data = np.array(raw)
x, y = data[:, 0], data[:, 1]
xrng = np.arange(-1.5, 1.5+0.05, 0.01)
yrng = np.arange(-1.5, 1.5+0.05, 0.01)
xx, yy = np.meshgrid( xrng, yrng )
zz = np.zeros(xx.shape)


def find_closest_index(val, arr):
    minidx = -1
    minval = 1000000
    for i in range(len(arr)):
        curval = np.abs( val - arr[i] )
        if curval < minval:
            minval = curval
            minidx = i
    return minidx

for i in range(len(x)):
    xix = find_closest_index( x[i], xrng )
    yix = find_closest_index( y[i], yrng )
    zz[xix,yix] += 1

np.save('intermediate/fig_2_xx.npy', xx)
np.save('intermediate/fig_2_yy.npy', yy)
np.save('intermediate/fig_2_zz.npy', zz)

# subfigure a : 3D representation

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = np.zeros(zz.shape)
p[:, :] = 1.1
ax.plot_surface(xx, yy, gaussian_filter(zz, 20), 
	cmap=plt.cm.coolwarm, linewidth=0)
ax.plot_surface(xx, yy, p, linewidth=0, 
	color='k', alpha=0.15)
plt.tick_params(labelsize=16)
ax.set_xlim(-1.50, 1.50)
ax.set_ylim(-1.50, 1.50)
ax.view_init(elev=35)
plt.savefig('figure_2_3d_repr.pdf')
# plt.show()

# subfigure b : contour 1

fig, ax = plt.subplots()
tmp = gaussian_filter(zz, 20)
tmp[tmp < 1.5]= 0
fig = plt.contour(tmp, colors=('#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E'))
plt.tick_params(labelsize=16)
ax.set_xticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
ax.set_yticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
cbar = plt.colorbar(fig)
cbar.ax.tick_params(labelsize=16)
plt.savefig('figure_2_contour_1.pdf')
# plt.show()

# subfigure c : contour 2

fig, ax = plt.subplots()
tmp = gaussian_filter(zz, 20)
tmp[tmp < 0.5]= 0
fig = plt.contour(tmp, colors=('#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E'))
plt.tick_params(labelsize=16)
ax.set_xticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
ax.set_yticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
cbar = plt.colorbar(fig)
cbar.ax.tick_params(labelsize=16)
plt.savefig('figure_2_contour_2.pdf')
# plt.show()

# save the figure to R so we can 
# get a persistence diagram

# import struct
# outputfile = 'intermediate/fig_2_data.bin'
# mat = data
# binfile = file(outputfile, 'wb')
# header = struct.pack('2I', mat.shape[0], mat.shape[1])
# header = struct.pack('2I', mat.shape[0], mat.shape[1])
# binfile.write(header)
# for i in range(mat.shape[1]):
#     data = struct.pack('%id' % mat.shape[0], *mat[:,i])
#     binfile.write(data)
# binfile.close()
