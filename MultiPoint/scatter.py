
# coding: utf-8

# # import modules

# In[ ]:


import numpy as np
from scipy import optimize
from numpy.linalg import eig

import MultiPoint as mp
from CORRcalc import *
from GAMcalc import *
from FEcalc import *

import propagator 
import wignerD as wd


# # define scattering functions

# In[ ]:


def scatter_pol(KV, N):
    FA=1.0  # leave this to 1.0 for homopolymers
    
    pset=propagator.prop_set(nlam=1)
    wigset = wd.wigner_d_vals()
    
    s2 = np.zeros((len(KV)),dtype=type(1+1j))
    for ind, K in enumerate(KV):
        s2[ind] = s2wlc(pset, N, FA, K, sequence=[0,0])/(N**2)
        
    return s2.real


# In[ ]:


def scatter_copol(KV, N, FA, CHI):
    pset=propagator.prop_set(nlam=1)

    gam2 = np.zeros((len(KV), 1),dtype=type(1+1j))
    for i, K in enumerate(KV):
        gam2[i] = gamma2(pset, N, FA, K, CHI)
        
    return 1./gam2.real


# In[ ]:


def scatter_sol(KV, N, FA, PHIP, CHIAB, CHIAS, CHIBS):
    pset=propagator.prop_set(nlam=1)
    mineig = np.zeros((len(KV), 1),dtype=type(1+1j))
    
    for i, K in enumerate(KV):
        mineig[i], _ = gamma2sol_smallereig(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS)
        
    return mineig.real

