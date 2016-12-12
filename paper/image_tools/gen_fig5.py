# -- Figure 5 : compilation of different tests --

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import seaborn as sns
from scipy.ndimage.filters import gaussian_filter

sns.set_style('whitegrid')

# e) kernel function

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.numpy2ri import numpy2ri

def rds_to_np(Rfile):
    ''' Convert .RData to be able to load
        the array into Python.
        Rfile := (str) location of file to
        translate to Python.
    '''
    raw = robjects.r['readRDS'](Rfile)
    return raw

# data = rds_to_np('intermediate/fig_5_kernel_pd.rds')
# data = np.array(data)
# data_feat_1 = data[:, 1:]
# x = data_feat_1[:, 0]
# y = data_feat_1[:, 1]

# xrng = np.arange(0.5, 1.7, 0.01)
# yrng = np.arange(0.5, 1.7, 0.01)
# xx, yy = np.meshgrid(xrng, yrng)
# zz = np.zeros(xx.shape)

# def find_closest_index(val, arr):
#     minidx = -1
#     minval = 1000000
#     for i in range(len(arr)):
#         curval = np.abs( val - arr[i] )
#         if curval < minval:
#             minval = curval
#             minidx = i

#     return minidx

# for i in range(len(x)):
#     xix = find_closest_index( x[i], xrng )
#     yix = find_closest_index( y[i], yrng )
#     zz[xix,yix] += 1

# plt.figure()
# tmp = gaussian_filter(zz.T, 1.5)
# plt.contour(xrng, yrng, tmp, linewidths=1, colors='k')
# fig = plt.contourf(xrng, yrng, tmp, cmap=plt.cm.jet)
# plt.tick_params(labelsize=28)
# plt.xlabel('Birth', fontsize=28)
# plt.ylabel('Death', fontsize=28)
# cbar = plt.colorbar(fig)
# cbar.ax.tick_params(labelsize=28)
# plt.tight_layout()
# plt.savefig('figure_5_kernel.pdf')
# # plt.show()

# # f) correlation function
# corr_fun = np.load('intermediate/fig_5_corr_func.npy')

# plt.figure()
# plt.plot(corr_fun, '-', linewidth=2.5)
# plt.tick_params(labelsize=28)
# plt.xlabel('Sequence', fontsize=28)
# plt.ylabel('CORR', fontsize=28)
# plt.tight_layout()
# plt.savefig('figure_5_corr_fun.pdf')
# # plt.show()

# # g) weighted kernel
# data = np.load('intermediate/fig_5_intensity.npy')
# plt.figure()
# fig = plt.pcolor(data.T, cmap=plt.cm.jet)
# plt.xlabel('Birth', fontsize=28)
# plt.ylabel('Death', fontsize=28)
# cbar = plt.colorbar(fig)
# cbar.ax.tick_params(labelsize=28)
# plt.tick_params(axis='both', which='major', labelsize=28)
# plt.tick_params(axis='both', which='minor', labelsize=28)
# plt.tight_layout()
# plt.savefig('figure_5_intensity_fun.pdf')
# # plt.show()

# h) persistent image
plt.figure()
data = np.load('intermediate/fig_5_pimage.npy')
data = data.reshape(13, 13)
fig = plt.pcolor(data, cmap=plt.cm.RdBu)
plt.xlabel('Birth', fontsize=28)
plt.ylabel('Persistence', fontsize=28)
cbar = plt.colorbar(fig)
cbar.ax.tick_params(labelsize=28)
plt.gca().set_xlim((0,13))
plt.gca().set_ylim((0,13))
plt.tick_params(axis='both', which='major', labelsize=28)
plt.tick_params(axis='both', which='minor', labelsize=28)
plt.tight_layout()
plt.savefig('figure_5_pimage_fun.pdf')
# plt.show()
