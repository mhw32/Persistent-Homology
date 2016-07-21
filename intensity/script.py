# File to run the intensity script.

# part 1 : Voronoi tests

from intensity import intensity_voronoi_test_suite, pimage_voronoi_test_suite
import numpy as np
import argparse

def run(mini=1, maxi=101):
    print("running instance with min %d, max %d" % (mini, maxi))
    intensity_stats_norm   = []
    intensity_stats_unnorm = []
    pimage_stats_norm      = []
    pimage_stats_unnorm    = []

    for iters in range(mini, maxi):
        base_file = 'data/voronoi/baseline%d.npy' % iters
        foam_file = 'data/voronoi/foam%d.npy' % iters

        print("iteration %d: unstandardized intensity calculation..." % iters)
        intensity_stats_unnorm.append(intensity_voronoi_test_suite(base_file, foam_file, False))
        print("iteration %d: standardized intensity calculation..." % iters)
        intensity_stats_norm.append(intensity_voronoi_test_suite(base_file, foam_file, True))

        print("iteration %d: unstandardized persistent image calculation..." % iters)
        pimage_stats_unnorm.append(pimage_voronoi_test_suite(base_file, foam_file, False))
        print("iteration %d: standardized persistent image calculation..." % iters)
        pimage_stats_norm.append(pimage_voronoi_test_suite(base_file, foam_file, True))

    np.save('output/voronoi/intensity_stats_unnorm.npy', intensity_stats_unnorm)
    np.save('output/voronoi/intensity_stats_norm.npy', intensity_stats_norm)
    np.save('output/voronoi/pimage_stats_unnorm.npy', pimage_stats_unnorm)
    np.save('output/voronoi/pimage_stats_norm.npy', pimage_stats_norm)

if __name__ == '__main__':
    # Construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('--min', type=int, default=1, help='minimum value')
    ap.add_argument('--max', type=int, default=101, help='maximum value')

    # parse args and run stuff
    args = vars(ap.parse_args())
    run(args['min'], args['max'])


