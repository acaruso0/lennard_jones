from scipy.constants import physical_constants


amu = physical_constants["atomic mass constant"]
k = physical_constants["Boltzmann constant"]
na = physical_constants["Avogadro constant"]

periodic_table = {
        "Ar": (18, 39.95) 
        }

#lj_interaction_table = {("Ar", "Ar"): (0.2381, 3.405e-10)}
lj_interaction_table = {("Ar", "Ar"): (1.65e-21, 3.4e-10)}
#lj_interaction_table = {("Ar", "Ar"): (0.185, 3.5)}
#lj_interaction_table = {("Ar", "Ar"): (0.294, 3.73)} # FAKE ARGON

