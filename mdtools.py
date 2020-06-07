import numpy as np
from scipy.stats import maxwell
from constants import k


def generate_initial_velocities(atoms, m, T):
    v = []
    for n in range(len(atoms)):
        a = np.sqrt((k[0]*T)/m[n])
        transl_vec = np.full(3, maxwell(scale=a).stats('m')/2)
        distr = maxwell.rvs(scale=a, size=3)
        v.append((distr - transl_vec)*2)
    return np.array(v)

def pbc_ortho(coords_copy, box, point=np.array([0,0,0])):
    coords_copy -= point
    coords_copy += np.divide(box, 2.)
    coords_copy -= (coords_copy // box)*box
    coords_copy -= np.divide(box, 2.)
    return coords_copy

