# File to run the intensity script.

# part 1 : Voronoi tests

from intensity import intensity_voronoi_test_suite, pimage_voronoi_test_suite

for iter in range(0, 100):
    base_file = 'data/voronoi/baseline%d.npy' % iter
    foam_file = 'data/voronoi/foam%d.npy' % iter

    intensity_voronoi_test_suite(base_file, foam_file, False)
    intensity_voronoi_test_suite(base_file, foam_file, True)

    pimage_voronoi_test_suite(base_file, foam_file, False)
    pimage_voronoi_test_suite(base_file, foam_file, True)
