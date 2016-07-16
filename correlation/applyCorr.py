from treeCorr import get_corr_stat, get_corr_func
from translate import rds_to_np
from translate import read_pure_foam, read_pure_baseline
from tools import normFoam, normVec
import numpy as np

# ---------------------------------------------------------------------------
# This returns the correlation statistic (mean of the absolute value of 
# correlation function) for both the voronoi foam data set
# and the simulation data set.

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

	# storage for this stuff
	base_stats = np.zeros(num_samples)
	foam_stats = np.zeros((num_percfil, num_samples))

	for i in range(num_samples):
		cur_base = base_data[i].T
		base_stats[i] = get_corr_stat(cur_base, L=35, Lc=-5, ngal=5000)
		
		for p in range(num_percfil):
			cur_foam = foam_data[p][i].T
			# get the stats store them.
			foam_stats[p, i] = get_corr_stat(cur_foam, L=35, Lc=-5, ngal=5000)

	return base_stats, foam_stats

def corr_simu_test_suite(cdm_file, wdm_file, normalize=False):
	# read all the data
	cdm_raw  = rds_to_np(cdm_file)
	if cdm_file == 'pure_simu/cdm_1.rds':
		cdm_data = [np.array( cdm_raw ).T]
	else:
		cdm_data = [np.array( i ).T for i in cdm_raw]

	if normalize:
		cdm_data = normVec(cdm_data)

	wdm_raw  = rds_to_np(wdm_file)
	if wdm_file == 'pure_simu/wdm_1.rds':
		wdm_data = [np.array( wdm_raw ).T]
	else:
		wdm_data = [np.array( i ).T for i in wdm_raw]

	if normalize:
		wdm_data = normVec(wdm_data)

	# define some useful constants
	num_samples = len(cdm_data)

	# storage for this stuff
	cdm_stats = np.zeros(num_samples)
	wdm_stats = np.zeros(num_samples)

	for i in range(num_samples):
		cur_cdm = cdm_data[i].T
		cdm_stats[i] = get_corr_stat(cur_cdm, L=100, Lc=0, ngal=cur_cdm.shape[0])
		
		cur_wdm = wdm_data[i].T
		wdm_stats[i] = get_corr_stat(cur_wdm, L=100, Lc=0, ngal=cur_wdm.shape[0])

	return cdm_stats, wdm_stats

# ---------------------------------------------------------------------------
# This returns the correlation function forboth the voronoi foam data set
# and the simulation data set.

def corr_voronoi_func_suite(base_file, foam_file, normalize=False):
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
	num_radius  = 100

	# storage for this stuff
	base_func = np.zeros((num_samples, num_radius))
	foam_func = np.zeros((num_percfil, num_samples, num_radius))

	for i in range(num_samples):
		cur_base = base_data[i].T
		tmp_func = get_corr_func(cur_base, L=35, Lc=-5, ngal=5000, max_r=num_radius)
		base_func[i, :tmp_func.shape[0]] = tmp_func
		
		for p in range(num_percfil):
			cur_foam = foam_data[p][i].T
			# get the stats store them.
			tmp_func = get_corr_func(cur_foam, L=35, Lc=-5, ngal=5000, max_r=num_radius)
			foam_func[p, i, :tmp_func.shape[0]] = tmp_func 

	return base_func, foam_func

def corr_simu_func_suite(cdm_file, wdm_file, normalize=False):
	# read all the data
	cdm_raw  = rds_to_np(cdm_file)
	if cdm_file == 'pure_simu/cdm_1.rds':
		cdm_data = [np.array( cdm_raw ).T]
	else:
		cdm_data = [np.array( i ).T for i in cdm_raw]

	if normalize:
		cdm_data = normVec(cdm_data)

	wdm_raw  = rds_to_np(wdm_file)
	if wdm_file == 'pure_simu/wdm_1.rds':
		wdm_data = [np.array( wdm_raw ).T]
	else:
		wdm_data = [np.array( i ).T for i in wdm_raw]

	if normalize:
		wdm_data = normVec(wdm_data)

	# define some useful constants
	num_samples = len(cdm_data)
	num_radius  = 100

	# storage for this stuff
	cdm_func = np.zeros((num_samples, num_radius))
	wdm_func = np.zeros((num_samples, num_radius))

	for i in range(num_samples):
		cur_cdm = cdm_data[i].T
		tmp_func = get_corr_func(cur_cdm, L=100, Lc=0, ngal=cur_cdm.shape[0], max_r=num_radius)
		cdm_func[i, :tmp_func.shape[0]] = tmp_func 
		
		cur_wdm = wdm_data[i].T
		tmp_func = get_corr_func(cur_wdm, L=100, Lc=0, ngal=cur_wdm.shape[0], max_r=num_radius)
		wdm_func[i, :tmp_func.shape[0]] = tmp_func 

	return cdm_func, wdm_func

if __name__ == '__main__':
	# for normalize in [False, True]:
	# 	all_base_corr = np.zeros((100, 15))
	# 	all_foam_corr = np.zeros((100, 9, 15))

	# 	for i in range(1, 101):
	# 		print('Operating on set %d' % i)
	# 		base_corr, foam_corr = corr_voronoi_test_suite( 'pure_data/baseline%d.rds' % i, 
	# 														'pure_data/foam%d.rds' % i,
	# 														normalize=normalize)
	# 		all_base_corr[i-1, :] = base_corr
	# 		all_foam_corr[i-1, :, :] = foam_corr

	# 	np.save('output/base_corr_norm(%d).npy' % normalize, all_base_corr)
	# 	np.save('output/foam_corr_norm(%d).npy' % normalize, all_foam_corr)

	# for normalize in [False, True]:
	# 	all_cdm_corr = np.zeros((4, 64))
	# 	all_wdm_corr = np.zeros((4, 64))

	# 	for i in np.arange(1,5)[::-1]:
	# 		print('Operating on set %d' % i)
	# 		cdm_corr, wdm_corr = corr_simu_test_suite(	'pure_simu/cdm_%d.rds' % i, 
	# 													'pure_simu/wdm_%d.rds' % i,
	# 													normalize=normalize)
	# 		all_cdm_corr[i-1, :i**3] = cdm_corr
	# 		all_wdm_corr[i-1, :i**3] = wdm_corr

	# 	np.save('output/cdm_corr_norm(%d).npy' % normalize, all_cdm_corr)
	# 	np.save('output/wdm_corr_norm(%d).npy' % normalize, all_wdm_corr)

	for normalize in [False, True]:
		all_base_corr = np.zeros((100, 15, 100))
		all_foam_corr = np.zeros((100, 9, 15, 100))

		for i in range(1, 101):
			print('Operating on set %d' % i)
			base_corr, foam_corr = corr_voronoi_func_suite( 'pure_data/baseline%d.rds' % i, 
															'pure_data/foam%d.rds' % i,
															normalize=normalize)
			all_base_corr[i-1, :, :] = base_corr
			all_foam_corr[i-1, :, :, :] = foam_corr

		np.save('output/funcs/base_corr_func_norm(%d).npy' % normalize, all_base_corr)
		np.save('output/funcs/foam_corr_func_norm(%d).npy' % normalize, all_foam_corr)

	for normalize in [False, True]:
		all_cdm_corr = np.zeros((4, 64, 100))
		all_wdm_corr = np.zeros((4, 64, 100))

		for i in np.arange(1,5)[::-1]:
			print('Operating on set %d' % i)
			cdm_corr, wdm_corr = corr_simu_func_suite(	'pure_simu/cdm_%d.rds' % i, 
														'pure_simu/wdm_%d.rds' % i,
														normalize=normalize)
			all_cdm_corr[i-1, :i**3, :] = cdm_corr
			all_wdm_corr[i-1, :i**3, :] = wdm_corr

		np.save('output/funcs/cdm_corr_func_norm(%d).npy' % normalize, all_cdm_corr)
		np.save('output/funcs/wdm_corr_func_norm(%d).npy' % normalize, all_wdm_corr)

