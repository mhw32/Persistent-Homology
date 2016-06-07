from treeCorr import get_corr_stat
from translate import rds_to_np
from translate import read_foam, read_baseline
import numpy as np

def corr_test_suite(base_file, foam_file):
	# read all the data
	base_raw = rds_to_np(base_file)
	base_data = read_baseline(base_raw)

	foam_raw = rds_to_np(foam_file)
	foam_data = read_foam(foam_raw)

	# define some useful constants
	num_samples = len(base_data)
	num_percfil = len(foam_data)
	num_dims    = len(np.unique(base_data[0][:, 0]))

	# storage for this stuff
	base_stats = np.zeros((num_samples, num_dims))
	foam_stats = np.zeros((num_percfil, num_samples, num_dims))

	for i in range(num_samples):
		cur_base = base_data[i]
		base_stats[i, :] = get_corr_stat(cur_base)
		
		for p in range(num_percfil):
			cur_foam = foam_data[p][i]
			# get the stats store them.
			foam_stats[p, i, :] = get_corr_stat(cur_foam)

	return base_stats, foam_stats

if __name__ == '__main__':
	all_base_corr = np.zeros((100, 15, 3))
	all_foam_corr = np.zeros((100, 9, 15, 3))

	for i in range(1, 101):
		print('Operating on set %d' % i)
		base_corr, foam_corr = corr_test_suite('data/baseline%d.rds' % i, 'data/foam%d.rds' % i)
		all_base_corr[i, :, :] = base_corr
		all_foam_corr[i, :, :, :] = foam_corr

	np.save('output/base_corr.npy', all_base_corr)
	np.save('output/foam_corr.npy', all_foam_corr)


