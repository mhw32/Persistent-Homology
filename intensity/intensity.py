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
    num = len(births)
    sum_intensity = 0
    
    for b, d in zip(births, deaths):
        intensity = (d - b) * (1 / tau**2) * gaussianKernel1D((x - b) / tau, h=np.std(births)) * gaussianKernel1D((y - d) / tau, h=np.std(deaths))
        # print(x, b, y, d, intensity)
        sum_intensity += intensity

    return sum_intensity

def intensityDiagFunc(diag, dim, delta, tau):
    diag = cleanDiag(diag)
    input = sliceDim(diag, dim)
    xvec, yvec = input[:, 1], input[:, 2]

    xlist = np.arange(min(xvec), max(xvec), delta)
    ylist = np.arange(min(yvec), max(yvec), delta)
    numx, numy = len(xlist), len(ylist)

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            if (y >= x):
                stats[i, j] = intensityEqn(x, y, xvec, yvec, tau)

    return (xlist, ylist, stats)
    # return stats

def intensityVecFunc(diagset, dim, delta, tau):
    f = lambda diag: intensityDiagFunc(diag, dim, delta, tau)
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
    elif t < b: return t**2 / float(b)
    else: return 1

def surfaceEqn(x, y, births, deaths):
    space = 0
    bias = np.max(deaths)
    for i, j in zip(births, deaths):
        space += weightFunc(j, b=bias) * gaussianKernel2D(x, y, mux=i, muy=j, hx=np.std(births), hy=np.std(deaths))
    return space      

def surfaceDiagFunc(diag, dim, delta):
    diag = cleanDiag(diag)
    diag = linearTransformDiag(diag)
    input = sliceDim(diag, dim)
    births, deaths = input[:, 1], input[:, 2]

    xlist = np.arange(min(births), max(births), delta)
    ylist = np.arange(min(deaths), max(deaths), delta)
    numx, numy = len(xlist), len(ylist)

    stats = np.zeros((numx, numy))
    for i in range(xlist.shape[0]):
        x = xlist[i]
        for j in range(ylist.shape[0]):
            y = ylist[j]
            stats[i, j] = surfaceEqn(x, y, births, deaths)

    return (xlist, ylist, stats)

# n, m should be divisble 
def surfaceToPixels(surf, n, m):
    # create image 
    rows, cols = surf.shape
    grid = np.zeros((n, m))
    grid_row = int(np.floor(rows / n))
    grid_col = int(np.floor(cols / n))

    # its persistence image is the collection of pixels I
    # is the double integral of the slice.
    for i in range(n):
        for j in range(m):
            grid[i,j] = np.sum(surf[grid_row*(i):grid_row*(i+1), 
                                    grid_col*(j):grid_col*(j+1)])

    return grid

def vectorize(A):
    return A.flatten()

def pimageDiagFunc(diag, delta, n, m):
    p0 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 0, 0.01), n, m))
    p1 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 1, 0.01), n, m))
    p2 = vectorize(surfaceToPixels(surfaceDiagFunc(diag, 2, 0.01), n, m))
    return np.concatenate([p0, p1, p2])

def pimageVecFunc(diagset, delta, n, m):
    f = lambda diag: pimageDiagFunc(diag, delta, n, m)
    return apply_to_vector(diagset, f)