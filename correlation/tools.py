import numpy as np

# Some tools.r converted to Python.

# clean a single persistence diagram
# def cleanDiag(diag):
# 	return diag[2:,:]

# # apply to each of the nsample ones
# def cleanVec(vec):
# 	return [cleanDiag(d) for d in vec]

# def cleanFoam(foam):
# 	cleaned = [cleanVec(v) for v in foam]
# 	return cleaned

# normalize a single persistence diagram
# def normDiag(diag):
# 	births, deaths = diag[:, 1], diag[:, 2]

# 	minbirths, mindeaths = min(births), min(deaths)
# 	maxbirths, maxdeaths = max(births), max(deaths)

# 	maxtotal = max([maxbirths, maxdeaths])
# 	mintotal = min([minbirths, mindeaths])

# 	newbirths = (births - mintotal) / float(maxtotal - mintotal)
# 	newdeaths = (deaths - mintotal) / float(maxtotal - mintotal)

# 	diag[:, 1], diag[:, 2] = newbirths, newdeaths
# 	return diag


def normSet(x):
	return (x - np.min(x)) / (np.max(x) - np.min(x))

# # normalize each of the nsample ones
def normVec(x):
	return [normSet(d) for d in x]

def normFoam(x):
	normed = [normVec(v) for v in x]
	return normed

