''' Translate R data to Python.  
'''

import rpy2.robjects as robjects
from rpy2.robjects.numpy2ri import numpy2ri
import numpy as np

def rds_to_np(Rfile):
	''' Convert .RData to be able to load 
		the array into Python.

		Rfile := (str) location of file to 
		translate to Python.
	''' 
	raw = robjects.r['readRDS'](Rfile)
	return raw

def read_baseline(raw):
	# convert raw matrices to numpy
	f = lambda x : np.array(x)
	return [f(x) for x in raw]

def read_foam(raw):
	return [read_baseline(x) for x in raw]

def read_fake_foam(raw):
	return [read_baseline(raw)]

def read_pure_baseline(raw):
	data = np.array(raw)
	new_data = np.zeros((data.shape[0] / 3, 3, data.shape[1]))
	for i in range(data.shape[0]):
		new_data[int(np.floor(i/3)),i%3,:] = data[i, :]
	return(list(new_data))

def read_pure_foam(raw):
	new_data = [] # can't be array, not all same length
	for i in range(len(raw)):
		new_data.append(read_pure_baseline(raw[i]))
	return new_data

def np_to_bin(inputfile, outputfile='/tmp/data.bin'):
	'''
	inputfile := npy file path
	outputfile := bin file path

	Read a numpy file, and write a simple 
	binary file containing two integers 'n' 
	and 'k' for rows and columns n times k floats 
	with the actual matrix which can be read 
	by any application or language that 
	can read binary.
	'''

	import struct
	import numpy as np

	# load from the file
	mat = np.load(inputfile)

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

def np_to_rds(inputfile, outputfile, robjname):
	''' Save numpy file to rds file.
	'''
	a = np.load(inputfile)
	a = np.array(a, dtype='float64')
	ro = numpy2ri(a)
	robjects.r.assign(robjname, ro)
	robjects.r("save(%s, file='%s', compress=TRUE)" % (robjname, outputfile))

