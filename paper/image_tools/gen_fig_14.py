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

NUMBER_TO_COLOR = 5

data = np.load('intermediate/fig_12_dists.npy')
data_dim = (data[-1, :, 0], data[-1, :, 1], data[-1, :, 2])

mat = np.zeros((3, NUMBER_TO_COLOR))

for i in range(3):
  d = np.argsort(data_dim[i])[:NUMBER_TO_COLOR]
  mat[i, :] = d

outputfile = "intermediate/fig_14_data.bin"
# create a binary file
binfile = file(outputfile, 'wb')
# and write out two integers with the row and column dimension
header = struct.pack('2I', mat.shape[0], mat.shape[1])
binfile.write(header)
# then loop over columns and write each
for i in range(mat.shape[1]):
    data = struct.pack('%id' % mat.shape[0], *mat[:,i])
    binfile.write(data)

binfile.close()
