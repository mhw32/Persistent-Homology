from helper import apply_to_vector
import numpy as np
from copy import copy
# import matplotlib.pyplot as plt

# clean a single persistence diagram
def cleanDiag(diag):
    return diag[2:,:]

# apply to each of the nsample ones
def cleanVec(vec):
    return [cleanDiag(d) for d in vec]

def cleanFoam(foam):
    cleaned = [cleanVec(v) for v in foam]
    return cleaned

def normSet(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))

def normDiag(x):
    x[:, 1] = normSet(x[:, 1])
    x[:, 2] = normSet(x[:, 2])
    return x

# # normalize each of the nsample ones
def normVec(x):
    return [normDiag(d) for d in x]

def normFoam(x):
    normed = [normVec(v) for v in x]
    return normed

# def plotDiag(diag):
#     plt.figure()
#     plt.plot(diag[:, 1], diag[:, 2], 'o', color='k')
#     plt.show()

def normalize1D(x):
    const = np.sum(x)
    return x / const

def map0to1(x):
    maxi = np.max(x)
    mini = np.min(x)
    return (x - mini) / maxi

# def plotIntensity(x, y, z):
#     (xgrid, ygrid), zgrid = np.meshgrid(x, y), z.T
#     zgrid_min, zgrid_max = zgrid.min(), zgrid.max()

#     plt.figure()
#     plt.pcolor(xgrid, ygrid, zgrid, cmap='jet', vmin=zgrid_min, vmax=zgrid_max)
#     # set the limits of the plot to the limits of the data
#     plt.axis([xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()])
#     plt.colorbar()    
#     plt.show()   

def safe_concatenate(x):
    y = np.concatenate(x) 
    if np.ndim(y) != 2:
        y = safe_concatenate(y)
    return y

def getDiagMinMax(base_data, foam_data, by_dim=True):
    data = np.concatenate(([np.concatenate(base_data), np.concatenate(np.concatenate(foam_data))]))
    num_dim = np.unique(base_data[0][:, 0]).shape[0]

    xmin, xmax = np.zeros(num_dim), np.zeros(num_dim)
    ymin, ymax = np.zeros(num_dim), np.zeros(num_dim)

    for d in range(num_dim):
        xmin[d] = np.min(sliceDim(data, d)[:, 1])
        xmax[d] = np.max(sliceDim(data, d)[:, 1])
        ymin[d] = np.min(sliceDim(data, d)[:, 2])
        ymax[d] = np.max(sliceDim(data, d)[:, 2])

    if not by_dim:
        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    return (xmin, xmax, ymin, ymax)        

def getDiagMinMaxSimu(cdm_data, wdm_data, by_dim=True):
    data = np.concatenate([np.concatenate(cdm_data), np.concatenate(wdm_data)])
    num_dim = np.unique(cdm_data[0][:, 0]).shape[0]

    xmin, xmax = np.zeros(num_dim), np.zeros(num_dim)
    ymin, ymax = np.zeros(num_dim), np.zeros(num_dim)

    for d in range(num_dim):
        xmin[d] = np.min(sliceDim(data, d)[:, 1])
        xmax[d] = np.max(sliceDim(data, d)[:, 1])
        ymin[d] = np.min(sliceDim(data, d)[:, 2])
        ymax[d] = np.max(sliceDim(data, d)[:, 2])

    if not by_dim:
        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    return (xmin, xmax, ymin, ymax)        


# -----------------------------------------------------------
# implementing a test statistic based on intensity
# images on persistence diagrams. 
# http://arxiv.org/pdf/1510.02502v1.pdf

def gaussianKernel1D(x, mu=0, h=1):
    k = (1 / (np.sqrt(2 * np.pi) * h)) * np.exp(-(x - mu)**2 / (2 * h**2))
    return k

def gaussianKernel2D(x, y, mux=0, muy=0, hx=1, hy=1):
    # (1 / (2 * np.pi * (hx**2 + hy**2))) 
    # i don't the constant and it doesn't really matter since all relative.
    k = np.exp(-((x - mux)**2 / (2 * hx**2) + (y - muy)**2 / (2 * hy**2)))
    return k

def gaussianKernel1DNoConst(x, mu=0, h=1):
    k = np.exp(-(x - mu)**2 / (h**2))
    return k

def sliceDim(diag, dim):
    return diag[diag[:, 0] == dim, :]

def intensityEqn(x, y, births, deaths, tau):
    intensity = (deaths - births) * (1 / tau**2) * gaussianKernel1D((x - births) / tau, h=np.std(births)) * gaussianKernel1D((y - deaths) / tau, h=np.std(deaths))
    return np.sum(intensity)

def intensityDiagFunc(diag, dim, gridnum, tau, xmin, xmax, ymin, ymax):
    inputs = sliceDim(diag, dim)
    xvec, yvec = inputs[:, 1], inputs[:, 2]

    xlist = np.linspace(xmin, xmax, gridnum)
    ylist = np.linspace(ymin, ymax, gridnum)
    numx, numy = len(xlist), len(ylist)

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            if (y >= x):
                stats[i, j] = intensityEqn(x, y, xvec, yvec, tau)

    return (xlist, ylist, stats)

def intensityVecDiagFunc(diag, dim, gridnum, tau, xmin, xmax, ymin, ymax):
    _, _, z = intensityDiagFunc(diag, dim, gridnum, tau, xmin, xmax, ymin, ymax)
    return vectorize(z)

def intensityVecFunc(diagset, dim, gridnum, tau, xmin, xmax, ymin, ymax):
    f = lambda diag: intensityVecDiagFunc(diag, dim, gridnum, tau, xmin, xmax, ymin, ymax)
    return apply_to_vector(diagset, f)
    
# -----------------------------------------------------------
# implementing a test statistic based on persistent images. 
# data --> persistence diagram --> surface --> image
# http://arxiv.org/pdf/1507.06217v3.pdf

def linearTransformDiag(diag):
     ''' rotate the diagram along the diagonal space '''
     newdiag = copy(diag)
     newdiag[:, 2] = newdiag[:, 2] - newdiag[:, 1]
     return newdiag

def weightFunc(s, b):
    t = copy(s)
    t[t <= 0] = 0
    t[t < b] = t[t < b] / float(b)
    t[t >= b] = 1
    return t

def surfaceEqn(x, y, births, deaths, tau, bias):
    space = weightFunc(deaths, b=bias) * (1 / tau**2) * gaussianKernel2D((x - births) / tau, (y - deaths) / tau, mux=0, muy=0, hx=np.std(births), hy=np.std(deaths))
    return np.sum(space)      

def surfaceDiagFunc(diag, dim, gridnum, tau, xmin, xmax, ymin, ymax):
    diag = linearTransformDiag(diag)
    inputs = sliceDim(diag, dim)
    births, deaths = inputs[:, 1], inputs[:, 2]

    xlist = np.linspace(xmin, xmax, gridnum)
    ylist = np.linspace(ymin, ymax, gridnum)
    numx, numy = len(xlist), len(ylist)

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            stats[i, j] = surfaceEqn(x, y, births, deaths, tau, bias=ymax)

    return (xlist, ylist, stats)

def safeGroup(totalsize, groupsize, loosefrac=0.5):
    ''' Given a number (assume a range up to it), how to 
        safely split into groups of size groupsize while
        losing as little information as possible.
    '''

    # divide and get the gractional part.
    integral = int(np.floor(totalsize / float(groupsize)))

    # breakpoints (naive)
    bkpts = range(0, totalsize, integral)

    # if there is any overflow, make it its own group if 
    # at least loosefrac * the size of the spacing
    if (totalsize - bkpts[-1]) >= (integral * loosefrac):
        bkpts.append(totalsize)

    return np.array(bkpts)

# n, m should be divisble 
def surfaceToPixels(surf, n):
    # create image 
    rows, cols = surf.shape

    if rows <= n or cols <= n:
        return surf

    rowgroupsizes = safeGroup(rows, n)
    colgroupsizes = safeGroup(cols, n)
    gridrows, gridcols = rowgroupsizes.shape[0] - 1, colgroupsizes.shape[0] - 1 
    grid = np.zeros((gridrows, gridcols))

    # its persistence image is the collection of pixels I
    # is the double integral of the slice.
    for i in range(gridrows):
        for j in range(gridcols):
            grid[i,j] = np.sum(surf[rowgroupsizes[i]:rowgroupsizes[i+1], colgroupsizes[j]:colgroupsizes[j+1]])

    return grid

def vectorize(A):
    return A.flatten()

def pimageDiagFunc(diag, gridnum, tau, n, xmin, xmax, ymin, ymax):
    p0 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 0, gridnum, tau, xmin, xmax, ymin, ymax)[2], n))
    p1 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 1, gridnum, tau, xmin, xmax, ymin, ymax)[2], n))
    p2 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 2, gridnum, tau, xmin, xmax, ymin, ymax)[2], n))
    return np.concatenate([p0, p1, p2])

def pimageVecFunc(diagset, gridnum, tau, n, xmin, xmax, ymin, ymax):
    f = lambda diag: pimageDiagFunc(diag, gridnum, tau, n, xmin, xmax, ymin, ymax)
    return apply_to_vector(diagset, f)

# -----------------------------------------------------------
# main scripts to apply the intensity tests.

def read_baseline(raw):
    return [np.array(x) for x in raw]

def read_foam(raw):
    return [read_baseline(x) for x in raw]

def read_pure_baseline(raw):
    data = np.array(raw)
    new_data = np.zeros((data.shape[0] / 3, 3, data.shape[1]))
    for i in range(data.shape[0]):
        new_data[int(np.floor(i/3)),i%3,:] = data[i, :]
    return(list(new_data))

def read_pure_foam(raw):
    new_data = [] # can't be array, not all same length
    for i in range(len(raw)):
        new_data.append(read_pure_baseline(raw[i]))
    return new_data

def process_voronoi_file(base_file, foam_file):
    # read all the data
    base_raw  = np.load(base_file)
    base_data = read_baseline(base_raw)
    foam_raw  = np.load(foam_file)
    foam_data = read_foam(foam_raw)

    return base_data, foam_data

def process_simu_file(cdm_file, wdm_file):
    # read all the data
    cdm_raw  = np.load(cdm_file)
    if type(cdm_raw[0]) == np.float:
        cdm_data = [np.array( cdm_raw )]
    else:
        cdm_data = [np.array( i ) for i in cdm_raw]

    wdm_raw  = np.load(wdm_file)
    if type(wdm_raw[0]) == np.float:
        wdm_data = [np.array( wdm_raw )]
    else:
        wdm_data = [np.array( i ) for i in wdm_raw]

    return cdm_data, wdm_data

def intensity_voronoi_test_suite(base_file, foam_file, normalize=False):
    # process the voronoi files into arrays
    base_data, foam_data = process_voronoi_file(base_file, foam_file)
    
    # remove known anomalies
    base_data = cleanVec(base_data)
    foam_data = cleanFoam(foam_data)

    # normalize the thing if need be
    if normalize:
        base_data = normVec(base_data)
        foam_data = normFoam(foam_data)

    xmin, xmax, ymin, ymax = getDiagMinMax(base_data, foam_data)

    num_samples, num_percfil, num_dim = len(base_data), len(foam_data), 3
    base_stats = np.zeros((num_dim, num_samples, 25**2))
    foam_stats = np.zeros((num_dim, num_percfil, num_samples, 25**2))

    for d in range(num_dim):
        # vectorized operations (do we need to transpose?)
        base_stats[d, :, :] = intensityVecFunc(base_data, d, 25, 0.1, xmin[d], xmax[d], ymin[d], ymax[d])

        # define some useful constants
        for p in range(num_percfil):
            foam_stats[d, p, :, :] = intensityVecFunc(foam_data[p], d, 25, 0.1, xmin[d], xmax[d], ymin[d], ymax[d])

    return base_stats, foam_stats

def intensity_simu_test_suite(cdm_file, wdm_file, dim, normalize=False):
    # process the simulation files into arrays
    cdm_data, wdm_data = process_simu_file(cdm_file, wdm_file)
    cdm_data, wdm_data = cleanVec(cdm_data), cleanVec(wdm_data)

    # normalize the thing if need be
    if normalize:
        cdm_data, wdm_data = normVec(cdm_data), normVec(wdm_data)

    xmin, xmax, ymin, ymax = getDiagMinMaxSimu(cdm_data, wdm_data)

    # define storage containers
    num_samples, num_dim = len(cdm_data), 3
    cdm_stats = np.zeros((num_dim, num_samples, 25**2))
    wdm_stats = np.zeros((num_dim, num_samples, 25**2))

    # calculate statistics
    for d in range(num_dim):
        cdm_stats[d, :, :] = intensityVecFunc(cdm_data, d, 25, 0.1, xmin[d], xmax[d], ymin[d], ymax[d])
        wdm_stats[d, :, :] = intensityVecFunc(wdm_data, d, 25, 0.1, xmin[d], xmax[d], ymin[d], ymax[d])

    return cdm_stats, wdm_stats

def pimage_voronoi_test_suite(base_file, foam_file, normalize=False):
    # process the voronoi files into data
    base_data, foam_data = process_voronoi_file(base_file, foam_file)    
    base_data = cleanVec(base_data)
    foam_data = cleanFoam(foam_data)

    # normalize the thing if need be
    if normalize:
        base_data = normVec(base_data)
        foam_data = normFoam(foam_data)

    # get important numbers for reshaping
    xmin, xmax, ymin, ymax = getDiagMinMax(base_data, foam_data, by_dim=False)

    # define some useful constants
    num_samples, num_percfil = len(base_data), len(foam_data)
    base_stats = pimageVecFunc(base_data, 25, 0.2, 10, xmin, xmax, ymin, ymax)
    foam_stats = np.zeros((num_percfil, num_samples, base_stats.shape[-1]))
    for p in range(num_percfil):
        foam_stats[p, :, :] = pimageVecFunc(foam_data[p], 25, 0.2, 10, xmin, xmax, ymin, ymax)

    return base_stats, foam_stats    

def pimage_simu_test_suite(cdm_file, wdm_file, normalize=False):
    # process the simulation files into arrays
    cdm_data, wdm_data = process_simu_file(cdm_file, wdm_file)
    cdm_data, wdm_data = cleanVec(cdm_data), cleanVec(wdm_data)
    
    # normalize the thing if need be
    if normalize:
        cdm_data, wdm_data = normVec(cdm_data), normVec(wdm_data)
    
    xmin, xmax, ymin, ymax = getDiagMinMaxSimu(cdm_data, wdm_data, by_dim=False)

    # storage for this stuff
    cdm_stats = pimageVecFunc(cdm_data, 25, 0.2, 10, xmin, xmax, ymin, ymax)
    wdm_stats = pimageVecFunc(wdm_data, 25, 0.2, 10, xmin, xmax, ymin, ymax)

    return cdm_stats, wdm_stats

# ----------------

def euclideanDist(X, Y):
    return np.sqrt(np.sum((X - Y)**2))

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def kernel_stat(X, Y, h):
    n, m = len(X), len(Y)
    nn_grid = cartesian((range(n), range(n)))
    nm_grid = cartesian((range(n), range(m)))
    mm_grid = cartesian((range(m), range(m)))

    # For each of the double for loops, loop through the grid.
    sum1 = np.sum(gaussianKernel1DNoConst(X[nn_grid[:, 0]], X[nn_grid[:, 1]], h))
    sum2 = np.sum(gaussianKernel1DNoConst(X[nm_grid[:, 0]], Y[nm_grid[:, 1]], h))
    sum3 = np.sum(gaussianKernel1DNoConst(Y[mm_grid[:, 0]], Y[mm_grid[:, 1]], h))

    # Now that we have all out sums, calculate the Tstat.
    T = float(1)/(n**2)*sum1 - float(2)/(m*n)*sum2 + float(1)/(m**2)*sum3    
    return T

def grid_search(fun, X, Y, mini=0.1, maxi=5.0, step=0.1):
    grid = np.arange(mini, maxi, step)
    maxp = mini
    maxval = np.abs(fun(X, Y, mini))
    for p in grid:
        curval = np.abs(fun(X, Y, p))
        if curval > maxval:
            maxval = curval
            maxp = p
    return maxp

def permute(X, Y):
    n = X.shape[0]
    Z = np.vstack((X, Y))
    np.random.shuffle(Z)
    return Z[:n], Z[n:]

# Permutation method:
# It estimates the p-value through sampling. Asymptotically correct.
# https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/
# N --> number of permutations to do (should be a big number).
def permutation_method(X, Y, N=10000, h=None):
    if h is None:
        h = grid_search(kernel_stat, X, Y, 0.1, 5.0, 0.1)
    
    p = 0 # estimated pval
    loss_orig = kernel_stat(X, Y, h=h)
    for _ in range(N):
        X, Y = permute(X, Y)
        loss_new = kernel_stat(X, Y, h=h)
        if (loss_new <= loss_orig):
            p += 1
    
    p = p / float(N)
    return p

# ------------

def reshape_data(data, by_dim=False):
    num_iters = len(data)
    if by_dim:
        num_dims = len(data[0][1])
        num_percfil = len(data[0][1][0])
        num_reps = len(data[0][1][0][0])
        num_size = len(data[0][1][0][0][0])
        base_data = np.zeros((num_iters, num_dims, num_reps, num_size))
        foam_data = np.zeros((num_iters, num_dims, num_percfil, num_reps, num_size))

        for i in range(num_iters):
            for j in range(num_dims):
                base_data[i, j, :, :] = data[i][0][j]
                foam_data[i, j, :, :, :] = data[i][1][j]  
    else:
        num_percfil = len(data[0][1])
        num_reps = len(data[0][1][0])
        num_size = len(data[0][1][0][0])
        base_data = np.zeros((num_iters, num_reps, num_size))
        foam_data = np.zeros((num_iters, num_percfil, num_reps, num_size))

        for i in range(num_iters):
            base_data[i, :, :] = data[i][0]
            foam_data[i, :, :, :] = data[i][1]

    return base_data, foam_data

def reshape_simu_data(data, by_dim=False):
    num_iters = 4
    if by_dim:
        num_dims = len(data[0][0])
        num_reps = 64 # hardcoded max
        num_size = len(data[0][0][0][0])
        cdm_data = np.zeros((num_iters, num_dims, num_reps, num_size))
        wdm_data = np.zeros((num_iters, num_dims, num_reps, num_size))

        for i in range(num_iters):
            cur_reps = len(data[i][0][0])
            cdm_data[i, :, :cur_reps, :] = data[i][0]
            wdm_data[i, :, :cur_reps, :] = data[i][1]
    else:
        num_reps = 64
        num_size = len(data[0][0][0])
        cdm_data = np.zeros((num_iters, num_reps, num_size))
        wdm_data = np.zeros((num_iters, num_reps, num_size))

        for i in range(num_iters):
            cur_reps = len(data[i][0])
            cdm_data[i, :cur_reps, :] = data[i][0]
            wdm_data[i, :cur_reps, :] = data[i][1]

    return cdm_data, wdm_data

# -- script to convert raw data into base/foam seps --
# import cPickle
# data_by_dim = ['intensity_stats_norm.pkl', 'intensity_stats_unnorm.pkl']
# data_no_dim = ['pimage_stats_norm.pkl', 'pimage_stats_unnorm.pkl']

# for dfile in data_by_dim:
#     data = cPickle.load(open(dfile, 'rb'))
#     base_data, foam_data = reshape_data(data, by_dim=True)
#     np.save(open('base_' + dfile[:-4] + '.npy', 'wb'), base_data)
#     np.save(open('foam_' + dfile[:-4] + '.npy', 'wb'), foam_data)

# for dfile in data_no_dim:
#     data = cPickle.load(open(dfile, 'rb'))
#     base_data, foam_data = reshape_data(data, by_dim=False)
#     np.save(open('base_' + dfile[:-4] + '.npy', 'wb'), base_data)
#     np.save(open('foam_' + dfile[:-4] + '.npy', 'wb'), foam_data)
# -- script to convert raw data into cdm/wdm seps --
# import cPickle
# data_by_dim = ['intensity_stats_norm.pkl', 'intensity_stats_unnorm.pkl']
# data_no_dim = ['pimage_stats_norm.pkl', 'pimage_stats_unnorm.pkl']

# for dfile in data_by_dim:
#     data = cPickle.load(open(dfile, 'rb'))
#     cdm_data, wdm_data = reshape_simu_data(data, by_dim=True)
#     np.save(open('cdm_' + dfile[:-4] + '.npy', 'wb'), cdm_data)
#     np.save(open('wdm_' + dfile[:-4] + '.npy', 'wb'), wdm_data)

# for dfile in data_no_dim:
#     data = cPickle.load(open(dfile, 'rb'))
#     cdm_data, wdm_data = reshape_simu_data(data, by_dim=False)
#     np.save(open('cdm_' + dfile[:-4] + '.npy', 'wb'), cdm_data)
#     np.save(open('wdm_' + dfile[:-4] + '.npy', 'wb'), wdm_data)
# -- end script --

# apply the permutation test onto the voronoi and the 
# simulation data. This wil return logged values. 

def voronoi_bydim_hypo_suite(base_stats, foam_stats):
    num_iters = base_stats.shape[0]
    num_dims  = base_stats.shape[1]
    num_percfil = foam_stats.shape[2]

    log_p_grid = np.zeros((num_iters, num_dims, num_percfil))
    for i in range(log_p_grid.shape[0]):
        print('iter: %d' % i)
        for d in range(num_dims):
            for j in range(num_percfil):
                log_p = np.log(permutation_method(base_stats[i, d, :, :], foam_stats[i, d, j, :, :], N=1000))
                log_p_grid[i, d, j] = log_p

    return log_p_grid

def simu_bydim_hypo_suite(cdm_stats, wdm_stats):
    num_iters = cdm_stats.shape[0]
    num_dims  = cdm_stats.shape[1]
    
    log_p_grid = np.zeros((num_iters, num_dims))
    for i in range(log_p_grid.shape[0]):
        print('iter: %d' % i)
        cur_length = (i+1)**3
        for d in range(num_dims):
            log_p = np.log(permutation_method(cdm_stats[i, d, :cur_length, :], wdm_stats[i, d, :cur_length, :], N=1000))
            log_p_grid[i, d] = log_p
    
    return log_p_grid

def voronoi_nodim_hypo_suite(base_stats, foam_stats):
    num_iters = base_stats.shape[0]
    num_percfil = foam_stats.shape[1]

    log_p_grid = np.zeros((num_iters, num_percfil))
    for i in range(log_p_grid.shape[0]):
        print('iter: %d' % i)
        for j in range(num_percfil):
            log_p = np.log(permutation_method(base_stats[i, :, :], foam_stats[i, j, :, :], N=1000))
            log_p_grid[i, j] = log_p

    return log_p_grid

def simu_nodim_hypo_suite(cdm_stats, wdm_stats):
    num_iters = cdm_stats.shape[0]
    log_p_grid = np.zeros(num_iters)
    for i in range(log_p_grid.shape[0]):
        print('iter: %d' % i)
        cur_length = (i+1)**3
        log_p = np.log(permutation_method(cdm_stats[i, :cur_length, :], wdm_stats[i, :cur_length, :], N=10000))
        log_p_grid[i] = log_p

    return log_p_grid