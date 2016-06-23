from treeCorr import get_corr_stat
from translate import rds_to_np
from translate import read_pure_foam, read_pure_baseline
from tools import normFoam, normVec
import numpy as np

def corr_voronoi_test_suite(base_file, foam_file, normalize=False):
	# read all the data
	base_raw  = rds_to_np(base_file)
	base_data = read_pure_baseline(base_raw)

	if normalize:
		base_data = normVec(base_data)

	foam_raw  = rds_to_np(foam_file)
	foam_data = read_pure_foam(foam_raw)
	
	if normalize:
		foam_data = normFoam(foam_data)

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

def corr_simu_test_suite(cdm_file, wdm_file, normalize=False):
	# read all the data
	cdm_raw  = rds_to_np(cdm_file)
	cdm_data = read_pure_baseline(cdm_raw)

	if normalize:
		cdm_data = normVec(cdm_data)

	wdm_raw  = rds_to_np(wdm_file)
	wdm_data = read_pure_baseline(wdm_raw)
	
	if normalize:
		wdm_data = normVec(wdm_data)

	# define some useful constants
	num_samples = len(cdm_data)
	num_dims    = len(np.unique(cdm_data[0][:, 0]))

	# storage for this stuff
	cdm_stats = np.zeros((num_samples, num_dims))
	wdm_stats = np.zeros((num_samples, num_dims))

	for i in range(num_samples):
		cur_cdm = cdm_data[i]
		cdm_stats[i, :] = get_corr_stat(cur_cdm)
		
		cur_wdm = wdm_data[i]
		wdm_stats[i, :] = get_corr_stat(cur_wdm)

	return cdm_stats, wdm_stats


if __name__ == '__main__':
	# for normalize in [True]:
	# 	all_base_corr = np.zeros((100, 15, 3))
	# 	all_foam_corr = np.zeros((100, 9, 15, 3))

	# 	for i in range(1, 101):
	# 		print('Operating on set %d' % i)
	# 		base_corr, foam_corr = corr_voronoi_test_suite( 'data/baseline%d.rds' % i, 
	# 														'data/foam%d.rds' % i,
	# 														normalize=normalize)
	# 		all_base_corr[i-1, :, :] = base_corr
	# 		all_foam_corr[i-1, :, :, :] = foam_corr

	# 	np.save('output/base_corr_norm(%d).npy' % normalize, all_base_corr)
	# 	np.save('output/foam_corr_norm(%d).npy' % normalize, all_foam_corr)

	for normalize in [False, True]:
		all_cdm_corr = np.zeros((4, 64, 3))
		all_wdm_corr = np.zeros((4, 64, 3))

		for i in range(1, 5):
			print('Operating on set %d' % i)
			cdm_corr, wdm_corr = corr_simu_test_suite(	'data_sim/cdm%d.rds' % i, 
														'data_sim/wdm%d.rds' % i,
														normalize=normalize)
			all_cdm_corr[i-1, :i**3, :] = cdm_corr
			all_wdm_corr[i-1, :i**3, :] = wdm_corr

		np.save('output/cdm_corr_norm(%d).npy' % normalize, all_cdm_corr)
		np.save('output/wdm_corr_norm(%d).npy' % normalize, all_wdm_corr)



