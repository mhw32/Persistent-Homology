from helper import apply_to_vector
import numpy as np
from copy import copy
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# implementing a test statistic based on intensity
# images on persistence diagrams. 
# http://arxiv.org/pdf/1510.02502v1.pdf

def gaussianKernel1D(x, mu=0, h=1):
    k = (1 / (np.sqrt(2 * np.pi) * h)) * np.exp(-(x - mu)**2 / (2 * h**2))
    return k

def gaussianKernel2D(x, y, mux=0, muy=0, h=1):
    k = (1 / (2 * np.pi * h)) * np.exp(-((x - mux)**2 * (y - muy)**2) / (2 * h**2))
    return k

def sliceDim(diag, dim):
    return diag[diag[:, 0] == dim, :]

def plotDiag(diag):
    plt.figure()
    plt.plot(diag[:, 1], diag[:, 2], 'o', color='k')
    plt.show()

def intensityEqn(x, y, births, deaths, tau, sigma):
    num = len(births)
    sum_intensity = 0
    
    for b, d in zip(births, deaths):
        intensity = (d - b) * (1 / tau**2) * gaussianKernel1D((x - b) / tau, h=sigma) * gaussianKernel1D((y - d) / tau, h=sigma)
        # print(x, b, y, d, intensity)
        sum_intensity += intensity

    return(sum_intensity)

def intensityDiagFunc(diag, dim, delta, tau, sigma):
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
                stats[i, j] = intensityEqn(x, y, xvec, yvec, tau, sigma)

    return (xlist, ylist, np.rot90(stats))

def intensityVecFunc(diagset, dim, delta, tau, sigma):
    f = lambda diag: intensityDiagFunc(diag, dim, delta, tau, sigma)
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

def makeSurface(diag, x, y):
    space = 0
    for i, j in zip(diag[:, 1:]):
        space += weightFunc(j, b=1) * gaussianKernel2D(x, y, mux=i, muy=j, h=1)
    return space      

def weightFunc(t, b):
    if t <= 0: return 0
    elif t < b: return t / float(b)
    else: return 1

def makePixels(surf, n):
    # create image 
    numrows, numcols = surf.shape
    img = np.zeros((numrows - n, numcols - n))

    # its persistence image is the collection of pixels I
    # is the double integral of the slice.
    for i in range(numrows - n):
        for j in range(numcols - n):
            img[i,j] = np.sum(surf[i:i+n, j:j+n])

    return img 