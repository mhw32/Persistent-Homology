# File to run the intensity script.

# part 1 : Voronoi tests

from intensity import intensity_voronoi_test_suite, pimage_voronoi_test_suite
import numpy as np

intensity_stats_norm   = []
intensity_stats_unnorm = []
pimage_stats_norm      = []
pimage_stats_unnorm    = []

for iter in range(1, 101):
    print("processing iteration %d" % iter)
    base_file = 'data/voronoi/baseline%d.npy' % iter
    foam_file = 'data/voronoi/foam%d.npy' % iter

    intensity_stats_unnorm.append(intensity_voronoi_test_suite(base_file, foam_file, False))
    intensity_stats_norm.append(intensity_voronoi_test_suite(base_file, foam_file, True))

    pimage_stats_unnorm.append(pimage_voronoi_test_suite(base_file, foam_file, False))
    pimage_stats_norm.append(pimage_voronoi_test_suite(base_file, foam_file, True))

np.save('output/voronoi/intensity_stats_unnorm.npy', intensity_stats_unnorm)
np.save('output/voronoi/intensity_stats_norm.npy', intensity_stats_norm)
np.save('output/voronoi/pimage_stats_unnorm.npy', pimage_stats_unnorm)
np.save('output/voronoi/pimage_stats_norm.npy', pimage_stats_norm)
