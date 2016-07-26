# File to run the intensity script.

# part 1 : Voronoi tests

from intensity import voronoi_bydim_hypo_suite, voronoi_nodim_hypo_suite
from intensity import simu_bydim_hypo_suite, simu_nodim_hypo_suite
import numpy as np
import cPickle 
import argparse

def run(testtype, norm):
    normstr = 'norm' if norm else 'unnorm' 

    if testtype == 'voronoi-bydim': 
        base_stats = np.load('output/voronoi/base_intensity_stats_%s.npy' % normstr)
        foam_stats = np.load('output/voronoi/foam_intensity_stats_%s.npy' % normstr)
        log_p_grid = voronoi_bydim_hypo_suite(base_stats, foam_stats)
    elif testtype == 'voronoi-nodim':
        base_stats = np.load('output/voronoi/base_pimage_stats_%s.npy' % normstr)
        foam_stats = np.load('output/voronoi/foam_pimage_stats_%s.npy' % normstr)
        log_p_grid = voronoi_nodim_hypo_suite(base_stats, foam_stats)
    elif testtype == 'simu-bydim':
        cdm_stats = np.load('output/simu/cdm_intensity_stats_%s.npy' % normstr)
        wdm_stats = np.load('output/simu/wdm_intensity_stats_%s.npy' % normstr)
        log_p_grid = simu_bydim_hypo_suite(cdm_stats, wdm_stats)
    elif testtype == 'simu-nodim':
        cdm_stats = np.load('output/simu/cdm_pimage_stats_%s.npy' % normstr)
        wdm_stats = np.load('output/simu/wdm_pimage_stats_%s.npy' % normstr)
        log_p_grid = simu_nodim_hypo_suite(cdm_stats, wdm_stats)
    else:
        return None

    return log_p_grid

if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('--type', type=str, default='voronoi-bydim', help='test type')
    ap.add_argument('--norm', type=bool, default=False, help='normalize test')

    # parse args and run stuff
    args = vars(ap.parse_args())
    run(args['type'], args['norm'])


