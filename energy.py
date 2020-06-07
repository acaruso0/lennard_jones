import mdtools
import numpy as np
from scipy import spatial
from constants import amu, k, lj_interaction_table


def lennard_jones(atom1, atoms, dists):
    norm = np.linalg.norm(dists, axis=1)
    params = np.array([lj_interaction_table[(atom1, atom2)] for atom2 in atoms])
    try:
        epsilon, sigma = params[:,0], params[:,1]
    except:
        return 0, np.array([0,0,0])
    ratio = np.divide(sigma, norm)
    energy = np.sum(4*epsilon*(ratio**12 - ratio**6))
    d_energy = (24*epsilon/norm)*(2*ratio**12 - ratio**6)
    gradient = np.sum(d_energy.reshape(-1,1)*(dists/norm.reshape(-1,1)), axis=0)
    return energy, gradient

def calculate_energy_grad(atoms, coords, box):
    energies, grads = [], []
    for n_coord, coord in enumerate(coords):
        coords_tmp = mdtools.pbc_ortho(coords.copy(), box, coord)
        kdtree = spatial.KDTree(coords_tmp)
        nbs = kdtree.query_ball_point(np.asarray([0,0,0]), 9e-10)
        nbs.remove(n_coord)
        energy, grad = lennard_jones(atoms[n_coord], atoms[nbs], coords_tmp[nbs])
        energies.append(energy)
        grads.append(grad)
    return np.array(energies), np.array(grads)

