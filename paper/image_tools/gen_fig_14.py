import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.numpy2ri import numpy2ri
import struct


def rds_to_np(Rfile):
  ''' 
  Convert .RData to be able to load 
  the array into Python.

  Rfile := (str) location of file to 
  translate to Python.
  ''' 
  raw = robjects.r['readRDS'](Rfile)
  return raw

NUMBER_TO_COLOR = 3

data = np.load('intermediate/fig_12_dists.npy')
data_dim = (data[-1, :, 0], data[-1, :, 1], data[-1, :, 2])

min_mat = np.zeros((3, NUMBER_TO_COLOR))
max_mat = np.zeros((3, NUMBER_TO_COLOR))

for i in range(3):
  d = np.argsort(data_dim[i])
  min_mat[i, :] = d[:NUMBER_TO_COLOR]
  max_mat[i, :] = d[-NUMBER_TO_COLOR:]

outputfile = "intermediate/fig_14_min_data_top_3.bin"
# create a binary file
binfile = file(outputfile, 'wb')
# and write out two integers with the row and column dimension
header = struct.pack('2I', min_mat.shape[0], min_mat.shape[1])
binfile.write(header)
# then loop over columns and write each
for i in range(min_mat.shape[1]):
    data = struct.pack('%id' % min_mat.shape[0], *min_mat[:,i])
    binfile.write(data)

binfile.close()

outputfile = "intermediate/fig_14_max_data_top_3.bin"
# create a binary file
binfile = file(outputfile, 'wb')
# and write out two integers with the row and column dimension
header = struct.pack('2I', max_mat.shape[0], max_mat.shape[1])
binfile.write(header)
# then loop over columns and write each
for i in range(max_mat.shape[1]):
    data = struct.pack('%id' % max_mat.shape[0], *max_mat[:,i])
    binfile.write(data)

binfile.close()
