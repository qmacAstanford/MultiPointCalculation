
# coding: utf-8

# ## Vertex functions
# 
# This code calculates the vertex functions from random-phase-approximation
# of copolymer melts.

# In[1]:

from CORRcalc import s2wlc, s2inverse, s3wlc, s4wlc, norm
from itertools import product
import propagator 
import numpy as np
from scipy import optimize


# In[2]:

get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
# import mpld3
# mpld3.enable_notebook()


# In[3]:

def r2(N):
    return N - 0.5*(1-np.exp(-2*N))


# In[4]:

def spinodal(pset, N, FA):
    CHI = 0
    K0 = 1/np.sqrt(r2(N))
    
    KS = optimize.fmin(lambda K: np.real(gamma2(pset, N, FA, K, CHI)), K0,                      disp=False)
    
    return KS


# ## Quadratic vertex

# In[ ]:

def gamma2(pset, N, FA, K, CHI):
    s2inv = s2inverse(pset, N, FA, K)

    D = [1,-1]    # sign indicator
    G = 0
    for I0, I1 in product([0,1], repeat=2):
        G += s2inv[I0, I1]*D[I0]*D[I1]
        
    return -2*CHI + N*G


# ## Cubic Vertex

# In[5]:

def gamma3(pset, N, FA, Ks):
    K1, K2, K3 = Ks
    if norm(K1+K2+K3) >= 1e-10:
        raise('Qs must add up to zero')
        
    if not (abs(norm(K1)-norm(K2)) < 1e-5         and abs(norm(K2)-norm(K3)) < 1e-5):
        raise('Qs must have same length')
    
    s3 = s3wlc(pset, N, FA, Ks)
    s2inv = s2inverse(pset, N, FA, norm(K1))
    
    val = 0
    for I0, I1, I2 in product([0,1], repeat=3):
        val -= s3[I0][I1][I2]* (s2inv[I0][0] - s2inv[I0][1])*                               (s2inv[I1][0] - s2inv[I1][1])*                               (s2inv[I2][0] - s2inv[I2][1])

    return val*(N**2)


# ## Quartic Vertex

# In[ ]:

def gamma4(pset, wigset, N, FA, Ks):
    K1, K2, K3, K4 = Ks
    if not (norm(K1) == norm(K2) == norm(K3) == norm(K4)):
        raise('Qs must have same length')
    
    K = norm(K1)
    K12 = norm(K1+K2)
    K13 = norm(K1+K3)
    K14 = norm(K1+K4)
    
    s4 = s4wlc(pset, wigset, N, FA, Ks)
    s31 = s3wlc(pset, N, FA, [K1, K2, -K1-K2])
    s32 = s3wlc(pset, N, FA, [K1, K3, -K1-K3])
    s33 = s3wlc(pset, N, FA, [K1, K4, -K1-K4])
    
    s2inv = s2inverse(pset, N, FA, K)
    s21inv = s2inverse(pset, N, FA, K12)
    s22inv = s2inverse(pset, N, FA, K13)
    s23inv = s2inverse(pset, N, FA, K14)
    
    G4 = np.zeros((2,2,2,2),dtype=type(1+1j))
    for a1, a2, a3, a4 in product([0,1], repeat=4):
        for I0, I1 in product([0,1], repeat=2):
            G4[a1][a2][a3][a4] +=                 s31[a1][a2][I0]*s31[a3][a4][I1]*s21inv[I0][I1] +                 s32[a1][a2][I0]*s32[a3][a4][I1]*s22inv[I0][I1] +                 s33[a1][a2][I0]*s33[a3][a4][I1]*s23inv[I0][I1]
    G4 -= s4
    
    val = 0
    for I0, I1, I2, I3 in product([0,1], repeat=4):
        val += G4[I0][I1][I2][I3] *                (s2inv[I0][0] - s2inv[I0][1])*                (s2inv[I1][0] - s2inv[I1][1])*                (s2inv[I2][0] - s2inv[I2][1])*                (s2inv[I3][0] - s2inv[I3][1])
                
    return val*(N**3)

