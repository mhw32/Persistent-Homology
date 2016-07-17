import numpy as np

def apply_to_foam(foam, func):
    return [apply_to_vector(vec, func) for vec in foam]

def apply_to_vector(vector, func):
    return np.apply_along_axis(func, 1, vector)