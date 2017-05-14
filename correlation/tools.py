from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np


def normSet(x):
	return (x - np.min(x)) / (np.max(x) - np.min(x))


# normalize each of the nsample ones
def normVec(x):
	return [normSet(d) for d in x]


def normFoam(x):
	normed = [normVec(v) for v in x]
	return normed
