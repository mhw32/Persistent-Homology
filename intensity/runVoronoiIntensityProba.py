from intensity import (
    voronoi_bydim_hypo_suite,
    voronoi_nodim_hypo_suite,
)
import numpy as np
import cPickle
import argparse


def run(dataDir, outDir, testType='voronoi-bydim', norm=False):
    normstr = 'norm' if norm else 'unnorm'

    if testType == 'voronoi-bydim':
        base_stats = np.load('{}/base_intensity_stats_{}.npy'.format(dataDir, normstr))
        foam_stats = np.load('{}/foam_intensity_stats_{}.npy'.format(dataDir, normstr))
        log_p_grid = voronoi_bydim_hypo_suite(base_stats, foam_stats)
    elif testType == 'voronoi-nodim':
        base_stats = np.load('{}/base_pimage_stats_{}.npy'.format(dataDir, normstr))
        foam_stats = np.load('{}/foam_pimage_stats_{}.npy'.format(dataDir, normstr))
        log_p_grid = voronoi_nodim_hypo_suite(base_stats, foam_stats)
    else:
        raise Exception('test type {} not recognized'.format(testType))
    np.save(open('{}/voronoi_proba_{}_{}.npy'.format(
        outDir, testType, normstr), 'wb'), log_p_grid)


if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('dataDir', type=str, help='input folder')
    ap.add_argument('outDir', type=str, help='input folder')
    ap.add_argument('--testType', type=str, default='voronoi-bydim', help='test type')
    ap.add_argument('--norm', type=bool, default=False, help='normalize test')

    # parse args and run stuff
    args = vars(ap.parse_args())
    run(args['dataDir'], args['outDir'],
        args['testType'], args['norm'])
