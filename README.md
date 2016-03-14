# Persistent-Homology

An ongoing research project to use statistical analysis to extrapolate topological structure from simulations of the observable Universe under certain sets of assumptions. Specific focus lies on distinguishing conceptualized theories of dark matter. The research will be published in the Royal Astronomy Journal.

## How to Run
 
Instructions are provided for the following tasks:

- Generating voronoi tessellations.
- Running the testing library. 
- Applying the library to real simulation studies.

### Voronoi Approximations

In order to gauge the effectiveness of the hypothesis tests, voronoi approximations are necessary since there are no measures of correctness for the real dataset. In one simulation, I generate a set of voronoi tessellations, deemed a voronoi group. Voronoi tessellations each individually have customizable filament percentages, wall percentages, clutter percentages, etc. When generating a Voronoi group, all percentages are held constant excep the filament percentage, which is varied from `0.1` to `0.9`. This variation is to mimic the differences between the observable Universe under warm and cold dark matter assumptions. A voronoi group takes the following settings:

- `N` : The number of particles per Voronoi foam.
- `Boxlim` : The volume of the foam i.e. `c(0, 50)`.
- `res` : The resolution ("fineness") of the foam.
- `perturb` : The amount of noise to add.
- `groupN` : The number of tessellations per group.
- `baseline` : Which percent filament to treat as the baseline i.e. `0.1`. 
- `nameId` : The name to give this group; it will be saved as this name.

To run the generation process, see `main_tests/run_voronoi.r`. Given some settings, the function to call is:

    >> source('VoronoiFoam.r')
    >> voronoi_compilation(N=N, Boxlim=Boxlim, res=res, perturb=perturb, groupN=groupN, baseline=baseline, nameId=i)

### Hypothesis Testing Library

With the testing library, we want to create two sets of structures: the `foam` and the `base`. The `base` represents the baseline for which we compare each structure in the `foam` to. In our case, we consider an independently generated base of 0.1 percFil and an alternative base of 0.9 percFil. We would like the tests to be repeatable, so it is run 100 times. See `cluster_scripts/run_vtests.r` for more details.

    >> source('testlib.r')
    >> foam <- readRDS(paste('./saved_states/large_set/foam', i, '.rds', sep=''))
    >> base <- readRDS(paste('./saved_states/large_set/baseline', i, '-0.1.rds', sep=''))
    >> test_wrapper(foam, base, paste(i, '0.1baseNormFalse', sep='-'), FALSE)

### Real Data

There are two options when trying to analyze the full Eagle dataset. Either one can use the KDE test, for which only 1 sample is needed, or one can slice the data set into different cubes in order to bootstrap a large enough set to use the traditional `testlib.r`. 

To load the Eagle data:

    >> source("process_eagle.r")
    >> cdm <- load_CDM()
    >> wdm <- load_WDM()

To create the respective persistence diagrams:

    >> res <- 2
    >> boxlim <- c(0, 100) 
    >> cdm_diag <- gridDiag(cdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    >> wdm_diag <- gridDiag(wdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)

All that remains is to calculate the p-value:

    >> pval <- ks::kde.test(cdm_diag$diagram, wdm_diag$diagram)$pvalue

Alternatively, slicing is done with the following procedure. From here the same testing framework is applicable:

    >> cdm_slices <- slice_cube_robust(cdm, 2)
    >> wdm_slices <- slice_cube_robust(wdm, 2)
    >> cdm_diags <- persistify_set(cdm_slices)
    >> wdm_diags <- persistify_set(wdm_slices)

## Repository Folders

- `cluster_scripts` : Scripts designed to run the simulations and topological analysis on the Yale Grace cluster. This uses the LSF platform so all scripts are `*.bsub`. Each `*.bsub` file has a corresponding `*.r` file.
- `main_tests` : Important scripts for utilizing the hypothesis testing framework to compare simulations.
- `simple_tests` : Development scripts for analyzing properties of Voronoi tessellations and persistent homologies.
- `saved_states` : A temporary database for results and figures prior to processing.
- `simulations` : A sample of the Eagle simulations of the observable Universe under WDM and CDM assumptions.
- `sketch` : Contains images created by an illustrator-type program.
- `images` : Holds images generated for the paper writing process as well as small demonstrations of topological analysis. Also holds results from visualizing hypothesis test p-values.
- `web` : A static web page to more-or-less showcase the work.
- `writeup` : Contains many different Overleaf repositories including the proposal, the complete writeup, kernel-density analysis, etc.

### Repository R Scripts

- `distance.r` : Contains functions for different distance measurements (Bottleneck, Wasserstein, etc.)
- `distribution.r` : Contains functions for deriving a statistic from the distribution and contour hypothesis tests.
- `euler.r` : Implementation of the Euler characteristic, which in theory is topologically invariant. 
- `localtest.r` : Functions using the local KDE test.
- `multiassign.r` : Short cut for doing multiple variable assignment in a single line.
- `NHST.r` : Implementation of randomization-style null hypothesis significance tests for persistence diagrams.
- `process_eagle.r` : Implementation of cubic slicing for Eagle data.
- `summarize.r` : Functions for processing Euler, Silhouette, and Landscape hypothesis tests. 
- `testlib.r` : Almost wrapper library for all the individual functions. This is where the API should hook in. 
- `tools.r` : Random (, possibly) useful functions.
- `twosample.r` : Functions for permutations testing (Gaussian Kernel).
- `Voronoi3Dfct.r` : Implementation of Voronoi tessellation generations with spacing derived from KNN. 
- `VoronoiFoam.r` : Wrapper to test `Voronoi3Dfct.r`. 

### Repository Python Scripts

- `parse.py` : Given the messy output of `testlib.r`, this is a robust way of organizing the data in a format that is easy to work with. For example:

        >> import parse
        >> start = '../saved_states/third_push/results-' 
        >> end   = '-0.1base.txt'
        >> filepaths = [start+str(i)+end for i in range(1, 21)]
        >> resArr = np.array([parse.parse(f) for f in filepaths])

- `parse.py` also contains two important functions: `parse.prepare1d` and `parse.prepare2d` that take `resArr` and further parse it into a plot-representable form. Certain tests are performed per dimension and require `prepare2d`; others only need `prepare1d`. 

        >> singles = ['all-silh', 'euler']
        >> doubles = ['indiv_silh', 'contour', 'global-kde']
        >> bighash = {}
        >> for characteristic in singles:
        >>     bighash[characteristic] = parse.prepare1d(resArr, characteristic)
        >> for characteristic in doubles:
        >>     for dim in [0,1,2]:
        >>         bighash[characteristic+'-dim-'+str(dim)] = parse.prepare2d(resArr, characteristic, dim)

## Contributing / Questions

Feel free to contribute to the code. If you have questions, email me at `me@mikewuis.me`. If you submit a pull request, I'd be happy to review it and most probably accept. Would love feedback too! If you spot a bug, let me know!

