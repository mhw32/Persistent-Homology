# File to run the intensity script.

# part 1 : Voronoi tests

from intensity import (
    intensity_simu_test_suite,
    pimage_simu_test_suite,
)
import numpy as np
import cPickle
import argparse


def run(dataDir, outDir, mini=1, maxi=5):
    print("running instance with min %d, max %d" % (mini, maxi))
    intensity_stats_norm   = []
    intensity_stats_unnorm = []
    pimage_stats_norm      = []
    pimage_stats_unnorm    = []

    for iters in range(mini, maxi):
        cdm_file = '{}/cdm_diags_{}.npy'.format(dataDir, iters)
        wdm_file = '{}/wdm_diags_{}.npy'.format(dataDir, iters)

        print("iteration {}: unstd intensity".format(iters))
        intensity_stats_unnorm.append(
            intensity_simu_test_suite(
                cdm_file,
                wdm_file,
                False,
            ),
        )
        print("iteration {}: std intensity".format(iters))
        intensity_stats_norm.append(
            intensity_simu_test_suite(
                cdm_file,
                wdm_file,
                True,
            ),
        )
        print("iteration {}: unstd persistent image".format(iters))
        pimage_stats_unnorm.append(
            pimage_simu_test_suite(
                cdm_file,
                wdm_file,
                False,
            ),
        )
        print("iteration {}: std persistent image".format(iters))
        pimage_stats_norm.append(
            pimage_simu_test_suite(
                cdm_file,
                wdm_file,
                True,
            ),
        )

    cPickle.dump(
        intensity_stats_unnorm,
        open('{}/intensity_stats_unnorm.pkl'.format(outDir), 'wb')
    )
    cPickle.dump(
        intensity_stats_norm,
        open('{}/intensity_stats_norm.pkl'.format(outDir), 'wb'),
    )
    cPickle.dump(
        pimage_stats_unnorm,
        open('{}/pimage_stats_unnorm.pkl'.format(outDir), 'wb'),
    )
    cPickle.dump(
        pimage_stats_norm,
        open('{}/pimage_stats_norm.pkl'.formt(outDir), 'wb'),
    )


if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('dataDir', type=str, help='input folder')
    ap.add_argument('outDir', type=str, help='input folder')
    ap.add_argument('--min', type=int, default=1, help='minimum value')
    ap.add_argument('--max', type=int, default=5, help='maximum value')

    # parse args and run stuff
    args = vars(ap.parse_args())
    run(args['dataDir'], args['outDir'],
        args['min'], args['max'])


