import pandas as pd
import numpy as np # NOT needed after line 6 is gone


def readxyz(filename):
    xyzfile = pd.read_csv(filename, sep='\s+', skiprows=2, header=None)
    atoms = xyzfile.iloc[:,0].to_numpy()
    coords = xyzfile.iloc[:,1:4].to_numpy()
    #return atoms, coords
    return np.full(10, 'Ar'), (np.random.rand(10, 3) - np.array([.5, .5, .5]))*2*10

def writexyz(atoms, coords):
    at_df = pd.DataFrame(atoms)
    coord_df = pd.DataFrame(coords)
    merged = pd.concat([at_df, coord_df], axis=1)
    with open("traj.xyz", 'a+') as outfile:
        outfile.write(F"{len(atoms)}\n\n")
        merged.to_csv(outfile, sep=' ', float_format='%.8f', header=False, index=False)
    return None

