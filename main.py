import iotools
import mdtools
import energy
import numpy as np
from constants import amu, na, periodic_table
import configparser
import thermostat

settings = "settings.ini"
config = configparser.ConfigParser()
config.read(settings)
xyzfile = config["system"]["xyzfile"]
ensemble = config["system"]["ensemble"]
box = tuple(float(i)*1e-10 for i in config["system"]["box"].strip('[]').split())
Ti = float(config["system"]["starting temperature"])
Tf = float(config["system"]["temperature"])
dt = float(config["simulation"]["dt"])
max_it = int(config["simulation"]["max_it"])
it = 0
if ensemble == "NVT":
    ts = thermostat.Berendsen(Tf, dt)

atoms, coords = iotools.readxyz(xyzfile)
dt *= 1e-15
m = np.array([periodic_table[atom][1] for atom in atoms])*amu[0]
iotools.writexyz(atoms, coords)

v = mdtools.generate_initial_velocities(atoms, m, Ti)
_, grad = energy.calculate_energy_grad(atoms, coords, box)
a = - grad / m.reshape(-1,1)
while it < max_it:
    if it%1000 == 0:
        print(it)
    coords += v*dt + 0.5*a*(dt**2)
    mdtools.pbc_ortho(coords, box)
    _, grad = energy.calculate_energy_grad(atoms, coords, box)
    if it%100 ==0:
        print((np.sum(_)/4184)*na[0])
    a_new = - grad / m.reshape(-1,1)
    v += 0.5*(a_new + a)*dt
    if ensemble == "NVT":
        T = ts.measure(v, m)
#        if it%1000 == 0:
#            print(T)
        v = ts.rescale(v, T)
    a = a_new
    if it%100 ==0:
        iotools.writexyz(atoms, coords)
    it += 1
