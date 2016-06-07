from treeCorr import get_corr_stat

def corr_test_suite(base_file, foam_file):
	# read all the data
	base_raw = rds_to_np(base_file)
	base_data = read_baseline(base_raw)

	foam_raw = rds_to_np(foam_file)
	foam_data = read_foam(foam_raw)

	# define some useful constants
	num_samples = len(base_data)
	num_percfil = len(foam_data[0])
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
