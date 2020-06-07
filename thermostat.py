import numpy as np
from constants import k


class Thermostat():
    def __init__(self, Tf):
        self.Tf = Tf

    def rescale(self):
        raise NotImplementedError

    def measure(self, v, m):
        v_2 = np.asarray([np.dot(vec, vec) for vec in v])
        KE = np.sum(m*v_2)
        return KE/k[0]*len(v)


class Scaling(Thermostat):
    def __init__(self, Tf):
        super().__init__(Tf)

    def rescale(self, v, T, Tf):
        self.c = np.sqrt(self.Tf/T)
        v *= self.c
        return v


class Berendsen(Thermostat):
    def __init__(self, Tf, dt, tau=0.5):
        super().__init__(Tf)
        self.dt = dt
        self.tau = tau

    def rescale(self, v, T, Tf):
        self.c = np.sqrt(1 - (self.dt/self.tau)*(1-self.T_f/T))
        v += self.c
        return v
    
