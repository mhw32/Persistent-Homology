# -- Figure 2 : toy model demo -- 

import numpy as np
from scipy.ndimage.filters import gaussian_filter
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.numpy2ri import numpy2ri
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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


raw_xy = rds_to_np('intermediate/fig_2_grid.rds')
raw_z = rds_to_np('intermediate/fig_2_dtm.rds')

data_xy = np.array(raw_xy).T
data_x = data_xy[:, 0]
data_y = data_xy[:, 1]
data_z = np.expand_dims(np.array(raw_z), axis=1)

size = int(np.sqrt(data_xy.shape[0]))
xx = data_x.reshape((size, size))
yy = data_y.reshape((size, size))
zz = data_z.reshape((size, size))

# subfigure a : 3D representation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = np.zeros(zz.shape)
p[:, :] = .8
ax.plot_surface(xx, yy, p, linewidth=0, 
    color='k', alpha=0.6)
ax.plot_surface(xx, yy, zz, 
    cmap=plt.cm.coolwarm, linewidth=0, alpha=0.7)
plt.tick_params(labelsize=16)
ax.set_xlim(-1.50, 1.50)
ax.set_ylim(-1.50, 1.50)
ax.view_init(elev=15, azim=35)
m = cm.ScalarMappable(cmap=cm.coolwarm)
m.set_array(zz)
cbar = plt.colorbar(m)
cbar.ax.tick_params(labelsize=16)
# plt.savefig('figure_2_3d_repr.pdf')
plt.show()

# subfigure b : contour 1
fig, ax = plt.subplots()
tmp = zz
tmp[tmp > 0.8]= 0
fig = plt.contour(tmp, colors=('#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E'))
plt.tick_params(labelsize=16)
ax.set_xticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
ax.set_yticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
cbar = plt.colorbar(fig)
cbar.ax.tick_params(labelsize=16)
cbar.ax.get_children()[0].set_linewidths([5]*7)
plt.savefig('figure_2_contour_1.pdf')
# plt.show()

# subfigure c : contour 2
fig, ax = plt.subplots()
tmp = zz
tmp[tmp > 0.2]= 0
fig = plt.contour(tmp, colors=('#44B3C2', '#F1A94E', '#E45641', '#5D4C46', '#7B8D8E'))
plt.tick_params(labelsize=16)
ax.set_xticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
ax.set_yticklabels([-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5])
cbar = plt.colorbar(fig)
cbar.ax.tick_params(labelsize=16)
cbar.ax.get_children()[0].set_linewidths([5]*7)
plt.savefig('figure_2_contour_2.pdf')
# plt.show()

