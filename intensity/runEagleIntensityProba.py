from intensity import (
    simu_bydim_hypo_suite,
    simu_nodim_hypo_suite,
)
import numpy as np
import cPickle
import argparse


def run(dataDir, outDir, testType='simu-bydim', norm=False):
    normstr = 'norm' if norm else 'unnorm'

    if testtype == 'simu-bydim':
        cdm_stats = np.load('{}/cdm_intensity_stats_{}.npy'.format(dataDir, normstr))
        wdm_stats = np.load('{}/wdm_intensity_stats_{}.npy'.format(dataDir, normstr))
        log_p_grid = simu_bydim_hypo_suite(cdm_stats, wdm_stats)
    elif testtype == 'simu-nodim':
        cdm_stats = np.load('{}/cdm_pimage_stats_{}.npy'.format(dataDir, normstr))
        wdm_stats = np.load('{}/wdm_pimage_stats_{}.npy'.format(dataDir, normstr))
        log_p_grid = simu_nodim_hypo_suite(cdm_stats, wdm_stats)
    else:
        raise Exception('test type {} not recognized'.format(testType))
    np.save(open('{}/simu_proba_{}_{}.npy'.format(
        outDir, testType, normstr), 'wb'), log_p_grid)


if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('dataDir', type=str, help='input folder')
    ap.add_argument('outDir', type=str, help='input folder')
    ap.add_argument('--testType', type=str, default='simu-bydim', help='test type')
    ap.add_argument('--norm', type=bool, default=False, help='normalize test')

    # parse args and run stuff
    args = vars(ap.parse_args())
    run(args['dataDir'], args['outDir'],
        args['testType'], args['norm'])
