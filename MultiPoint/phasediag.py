from .GAMcalc import *
from .FEcalc import *
from . import propagator 
from . import wignerD as wd

import numpy as np
from scipy import optimize
from numpy.linalg import eig

filename = 'data/phasediag'
f1 = open(filename, 'w')

#NV = [1e3, 1e2, 10, 1]
#FAV = np.linspace(.15, .499, 21)
#PHIPV = [0.99, 0.75, 0.5, 0.25]

NV = [1e3, 1]
FAV = [0.15, 0.499]
PHIPV = [0.99, 0.5]

for N in NV:
    for PHIP in PHIPV:
        CHIsolsV = np.zeros((len(FAV), 1))
        CHIsol13V = np.zeros((len(FAV), 1))
        CHIsol36V = np.zeros((len(FAV), 1))

        for i, FA in enumerate(FAV):
            CHIsolsV[i], CHIsol13V[i], CHIsol36V[i] = chiootsol(N, FA, PHIP)
            print(PHIP, N, FA, CHIsolsV[i], CHIsol13V[i], CHIsol36V[i])
            f1.write('%.2f, %.2f, %.4f, %.4f, %.4f, %.4f\n' %(PHIP, N, FA,
                                                              CHIsolsV[i]*N*PHIP,
                                                              CHIsol13V[i]*N*PHIP,
                                                              CHIsol36V[i]*N*PHIP))

PHIP = 1.0            
for N in NV:
    CHIsV = np.zeros((len(FAV), 1))
    CHI13V = np.zeros((len(FAV), 1))
    CHI36V = np.zeros((len(FAV), 1))
    for i, FA in enumerate(FAV):
        CHIsV[i] = chis(N, FA)
        CHI13V[i], CHI36V[i] = chioot(N, FA)
        print(N, FA, CHIsV[i], CHI13V[i], CHI36V[i])
        f1.write('%.2f, %.2f, %.4f, %.4f, %.4f, %.4f\n' %(PHIP, N, FA,
                                                          CHIsV[i]*N,
                                                          CHI13V[i]*N,
                                                          CHI36V[i]*N))

f1.close()
