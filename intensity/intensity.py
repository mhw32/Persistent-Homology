from helper import apply_to_vector
import numpy as np
from copy import copy
import matplotlib.pyplot as plt

# clean a single persistence diagram
def cleanDiag(diag):
    return diag[2:,:]

# apply to each of the nsample ones
def cleanVec(vec):
    return [cleanDiag(d) for d in vec]

def cleanFoam(foam):
    cleaned = [cleanVec(v) for v in foam]
    return cleaned

def plotDiag(diag):
    plt.figure()
    plt.plot(diag[:, 1], diag[:, 2], 'o', color='k')
    plt.show()

def normalize1D(x):
    const = np.sum(x)
    return x / const

def map0to1(x):
    maxi = np.max(x)
    mini = np.min(x)
    return (x - mini) / maxi

def plotIntensity(x, y, z):
    (xgrid, ygrid), zgrid = np.meshgrid(x, y), z.T
    zgrid_min, zgrid_max = zgrid.min(), zgrid.max()

    plt.figure()
    plt.pcolor(xgrid, ygrid, zgrid, cmap='jet', vmin=zgrid_min, vmax=zgrid_max)
    # set the limits of the plot to the limits of the data
    plt.axis([xgrid.min(), xgrid.max(), ygrid.min(), ygrid.max()])
    plt.colorbar()    
    plt.show()   

def safe_concatenate(x):
    y = np.concatenate(x) 
    if np.ndim(y) != 2:
        y = safe_concatenate(y)
    return y

def getDiagMinMax(base_data, foam_data, by_dim=True):
    data = np.vstack([np.concatenate(base_data), np.concatenate(np.concatenate(foam_data))])
    num_dim = np.unique(base_data[0][:, 0]).shape[0]

    xmin, xmax = np.zeros(num_dim), np.zeros(num_dim)
    ymin, ymax = np.zeros(num_dim), np.zeros(num_dim)

    for d in range(num_dim):
        xmin[d] = np.min(sliceDim(data, d), 1)
        xmax[d] = np.max(sliceDim(data, d), 1)
        ymin[d] = np.min(sliceDim(data, d), 2)
        ymax[d] = np.max(sliceDim(data, d), 2)

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

def sliceDim(diag, dim):
    return diag[diag[:, 0] == dim, :]

def intensityEqn(x, y, births, deaths, tau):
    sum_intensity = 0
    for b, d in zip(births, deaths):
        intensity = (d - b) * (1 / tau**2) * gaussianKernel1D((x - b) / tau, h=np.std(births)) * gaussianKernel1D((y - d) / tau, h=np.std(deaths))
        sum_intensity += intensity

    return sum_intensity

def intensityDiagFunc(diag, dim, gridnum, tau):
    diag = cleanDiag(diag)
    inputs = sliceDim(diag, dim)
    xvec, yvec = inputs[:, 1], inputs[:, 2]

    xlist = np.linspace(xvec.min(), xvec.max(), gridnum)
    ylist = np.linspace(yvec.min(), yvec.max(), gridnum)
    numx, numy = len(xlist), len(ylist)

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            if (y >= x):
                stats[i, j] = intensityEqn(x, y, xvec, yvec, tau)

    return (xlist, ylist, stats)

def intensityVecDiagFunc(diag, dim, gridnum, tau):
    _, _, z = intensityDiagFunc(diag, dim, gridnum, tau)
    return vectorize(z)

def intensityVecFunc(diagset, dim, gridnum, tau):
    f = lambda diag: intensityVecDiagFunc(diag, dim, gridnum, tau)
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

def weightFunc(t, b):
    if t <= 0: return 0
    elif t < b: return t / float(b)
    else: return 1

def surfaceEqn(x, y, births, deaths, tau, bias):
    space = 0
    for i, j in zip(births, deaths):
        space += weightFunc(j, b=bias) * (1 / tau**2) * gaussianKernel2D((x - i) / tau, (y - j) / tau, mux=0, muy=0, hx=np.std(births), hy=np.std(deaths))
    return space      

def surfaceDiagFunc(diag, dim, gridnum, tau, bias=None):
    diag = cleanDiag(diag)
    diag = linearTransformDiag(diag)
    inputs = sliceDim(diag, dim)
    births, deaths = inputs[:, 1], inputs[:, 2]

    xlist = np.linspace(births.min(), births.max(), gridnum)
    ylist = np.linspace(deaths.min(), deaths.max(), gridnum)
    numx, numy = len(xlist), len(ylist)
    if bias is None: bias = deaths.max()

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            stats[i, j] = surfaceEqn(x, y, births, deaths, tau, bias=bias)

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

def pimageDiagFunc(diag, gridnum, tau, bias, n):
    p0 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 0, gridnum, tau, bias)[2], n))
    p1 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 1, gridnum, tau, bias)[2], n))
    p2 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 2, gridnum, tau, bias)[2], n))
    return np.concatenate([p0, p1, p2])

def pimageVecFunc(diagset, gridnum, tau, bias, n):
    f = lambda diag: pimageDiagFunc(diag, gridnum, tau, bias, n)
    return apply_to_vector(diagset, f)

# -----------------------------------------------------------
# main scripts to apply the intensity tests.
import sys
sys.path.append('../correlation')
from translate import rds_to_np
from translate import read_pure_foam, read_pure_baseline
from translate import read_foam, read_baseline
from tools import normFoam, normVec

def process_voronoi_file(base_file, foam_file, normalize):
    # read all the data
    base_raw  = rds_to_np(base_file)
    base_data = read_baseline(base_raw)

    if normalize:
        base_data = normVec(base_data)

    foam_raw  = rds_to_np(foam_file)
    foam_data = read_foam(foam_raw)

    if normalize:
        foam_data = normFoam(foam_data)
    
    return base_data, foam_data

def process_simu_file(cdm_file, wdm_file, normalize):
    # read all the data
    cdm_raw  = rds_to_np(cdm_file)
    if type(cdm_raw[0]) == np.float:
        cdm_data = [np.array( cdm_raw ).T]
    else:
        cdm_data = [np.array( i ).T for i in cdm_raw]

    if normalize:
        cdm_data = normVec(cdm_data)

    wdm_raw  = rds_to_np(wdm_file)
    if type(wdm_raw[0]) == np.float:
        wdm_data = [np.array( wdm_raw ).T]
    else:
        wdm_data = [np.array( i ).T for i in wdm_raw]

    if normalize:
        wdm_data = normVec(wdm_data)
    
    return cdm_data, wdm_data

def intensity_voronoi_test_suite(base_file, foam_file, normalize=False):
    # process the voronoi files into arrays
    base_data, foam_data = process_voronoi_file(base_file, foam_file, normalize)
    # xmin, xmax, ymin, ymax = getDiagMinMax(base_data, foam_data)

    num_samples, num_percfil, num_dim = base_data.shape[0], foam_data.shape[0], base_data.shape[1] 
    base_stats = np.zeros((num_dim, num_samples, 30**2))
    foam_stats = np.zeros((num_dim, num_percfil, num_samples, 30**2))

    for d in range(num_dim):
        # vectorized operations (do we need to transpose?)
        base_stats[d, :, :] = intensityVecFunc(base_data, d, 30, 0.1)

        # define some useful constants
        for p in range(num_percfil):
            foam_stats[d, p, :, :] = intensityVecFunc(foam_data[p, :], d, 30, 0.1)

    return base_stats, foam_stats

def intensity_simu_test_suite(cdm_file, wdm_file, dim, normalize=False):
    # process the simulation files into arrays
    cdm_data, wdm_data = process_simu_file(cdm_file, wdm_file, normalize)
    # xmin, xmax, ymin, ymax = getDiagMinMax(cdm_data, wdm_data)

    # define storage containers
    num_samples, num_dim = cdm_data.shape[0], cdm_data.shape[1]
    cdm_stats = np.zeros((num_dim, num_samples, 30**2))
    wdm_stats = np.zeros((num_dim, num_samples, 30**2))

    # calculate statistics
    for d in range(num_dim):
        cdm_stats[d, :, :] = intensityVecFunc(cdm_data, d, 30, 0.1)
        wdm_stats[d, :, :] = intensityVecFunc(wdm_data, d, 30, 0.1)

    return cdm_stats, wdm_stats

def pimage_voronoi_test_suite(base_file, foam_file, normalize=False):
    # process the voronoi files into data
    base_data, foam_data = process_voronoi_file(base_file, foam_file, normalize)
    _, _, _, ymax = getDiagMinMax(base_data, foam_data, by_dim=False)

    # define some useful constants
    num_samples, num_percfil = len(base_data), len(foam_data)
    base_stats = pimageVecFunc(base_data, 30, 0.1, ymax, 10)
    foam_stats = np.zeros((num_percfil, num_samples, 30**2))
    for p in range(num_percfil):
        foam_stats[p, :, :] = pimageVecFunc(foam_data[p, :], 30, 0.1, ymax, 10)

    return base_stats, foam_stats    

def pimage_simu_test_suite(cdm_file, wdm_file, normalize=False):
    # process the simulation files into arrays
    cdm_data, wdm_data = process_simu_file(cdm_file, wdm_file, normalize)
    _, _, _, ymax = getDiagMinMax(cdm_data, wdm_data, by_dim=False)

    # storage for this stuff
    cdm_stats = pimageVecFunc(cdm_data, 30, 0.1, ymax, 30)
    wdm_stats = pimageVecFunc(wdm_data, 30, 0.1, ymax, 30)

    return cdm_stats, wdm_stats