
# coding: utf-8

# In[1]:


import numpy as np
from scipy import optimize
from GAMcalc import *

import propagator
import wignerD as wd
from numpy.linalg import eig


# In[2]:


def set3Ks(K):
    k1 = np.array([1,0,0])*K
    k2 = np.array([-0.5,0.5*np.sqrt(3),0])*K
    k3 = -k1-k2
    return [k1, k2, k3]

def set4Ks1(K):
    k1=np.array([1,0,0])*K
    k2=np.array([1,0,0])*K
    k3, k4 = -k1, -k2
    return [k1, k2, k3, k4]

def set4Ks2(K):
    k1=np.array([1,0,0])*K
    k2=np.array([0.5,0.5*np.sqrt(3),0])*K
    k3, k4 = -k1, -k2
    return [k1, k2, k3, k4]

def set4Ks3(K):
    k1=np.array([1,0,0])*K
    k2=np.array([0,1,0])*K
    k3, k4 = -k1, -k2
    return [k1, k2, k3, k4]

def set4Ks4(K):
    k1=np.array([-1,0,1])/np.sqrt(2)*K
    k2=np.array([-1,0,-1])/np.sqrt(2)*K
    k3=np.array([1,1,0])/np.sqrt(2)*K
    k4=np.array([1,-1,0])/np.sqrt(2)*K
    return [k1, k2, k3, k4]


# ### 0.4 Minimizing functions

# In[3]:


def chis(N, FA):
    pset=propagator.prop_set(nlam=1)
    Ks = spinodal(pset, N, FA)
    chis = gamma2(pset, N, FA, Ks, 0).real/2
    
    return chis


# In[4]:


def gamma2sol_biggereig(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS):
    gam2sol = gamma2sol(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS)
    eigv, eigvec = eig(gam2sol)
    return eigv.max().real, eigvec[:,eigv.argmax()].real


# In[5]:


def gamma2sol_smallereig(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS):
    gam2sol = gamma2sol(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS)
    eigv, eigvec = eig(gam2sol)
    return eigv.min().real, eigvec[:,eigv.argmin()].real


# In[6]:


def gamma2sol_mineig(pset, N, FA, PHIP, CHIAB, CHIAS, CHIBS):
    f = lambda K: gamma2sol_smallereig(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS)[0]
    
    # solve for critical wavemode
    Ks = optimize.fminbound(f, 1e-3/np.sqrt(N), 1e2/np.sqrt(r2(N)), disp=False)
    _, eigvmin = gamma2sol_smallereig(pset, N, FA, PHIP, Ks, CHIAB, CHIAS, CHIBS)
    
    return f(Ks), Ks, eigvmin


# In[7]:


def gamma2sol_spinodal(pset, N, FA, PHIP):
    # check this, this seems to effect the initialization of the below solving for eqn.
    CHIAS, CHIBS = 0, 0
    
    mineig = lambda CHIAB: gamma2sol_mineig(pset, N, FA, PHIP, CHIAB, CHIAS, CHIBS)[0]
    CHIABs = optimize.brentq(mineig, 1e-2/N/PHIP, 500/N/PHIP, xtol=1e-2/N/PHIP)
    
    # find Ks at spinodal
    _, Kss, eigvmin = gamma2sol_mineig(pset, N, FA, PHIP, CHIABs, CHIAS, CHIBS)
    return CHIABs, Kss, eigvmin


# ### Calculate free energies

# In[8]:


def FE1(A, CHI, gam41):
    pi = 3.14159
    a = -(2*pi)**(-3)*2*CHI
    b = 0
    c = (2*pi)**(-9)*(1/4)*gam41.real

    return a*A**2 + b*A**3 + c*A**4

def FE3(A, CHI, gam3, gam41, gam42):
    pi = 3.14159
    a = -(2*pi)**(-3)*2*CHI
    b = (2*pi)**(-6)*(2/3/np.sqrt(3))*gam3.real
    c = (2*pi)**(-9)*(1/12)*(gam41+4*gam42).real
    
    return a*A**2 + b*A**3 + c*A**4

def FE6(A, CHI, gam3, gam41, gam42, gam43, gam44):
    pi = 3.14159
    a = -(2*pi)**(-3)*2*CHI
    b = (2*pi)**(-6)*(4/3/np.sqrt(6))*gam3.real
    c = (2*pi)**(-9)*(1/24)*(gam41+8*gam42+2*gam43+4*gam44).real
    
    return a*A**2 + b*A**3 + c*A**4


# In[9]:


def chioot(N, FA):
    pset=propagator.prop_set(nlam=5)
    wigset = wd.wigner_d_vals()
    Ks = spinodal(pset, N, FA)
    
    # calculate vertices
    gam3 = gamma3(pset, N, FA, set3Ks(Ks)).real
    gam41 = gamma4(pset, wigset, N, FA, set4Ks1(Ks)).real
    gam42 = gamma4(pset, wigset, N, FA, set4Ks2(Ks)).real
    gam43 = gamma4(pset, wigset, N, FA, set4Ks3(Ks)).real
    gam44 = gamma4(pset, wigset, N, FA, set4Ks4(Ks)).real
    
    # calculate oot temperatures
    A0 = 1e-1
    FE13sq = lambda CHI:     (optimize.fmin(lambda A: FE1(A, CHI, gam41),
                  A0, disp=False, full_output=1)[1] - \
    optimize.fmin(lambda A: FE3(A, CHI, gam3, gam41, gam42),
                  A0, disp=False, full_output=1)[1])**2
    
    FE36sq = lambda CHI:     (optimize.fmin(lambda A: FE3(A, CHI, gam3, gam41, gam42),
                  A0, disp=False, full_output=1)[1] -\
    optimize.fmin(lambda A: FE6(A, CHI, gam3, gam41, gam42, gam43, gam44),
                  A0, disp=False, full_output=1)[1])**2

    chi13 = optimize.fminbound(FE13sq, 0/N, 100/N)
    chi36 = optimize.fminbound(FE36sq, 0/N, 100/N)
    
    return chi13, chi36


# ### Free energy functions for solutions

# In[10]:


def FE1sol(A1, A2, gam2sol, gam41sol):
    # note here the amplitudes are in directions [1, 0] and [0, 1]
    # i.e. psi = [A1, A2]*sin(q*x)
    pi = 3.14159
    FE = 0
    A = [A1, A2]
    for I1, I2 in product([0,1], repeat=2):
        FE += (2*pi)**(-3)*gam2sol[I1][I2]*A[I1]*A[I2]

    for I1, I2, I3, I4 in product([0,1], repeat=4):
        FE += (2*pi)**(-9)*(1/4)*        gam41sol[I1][I2][I3][I4]*A[I1]*A[I2]*A[I3]*A[I4].real
        
    return FE

def FE3sol(A1, A2, gam2sol, gam3sol, gam41sol, gam42sol):
    pi = 3.14159
    FE = 0
    A = [A1, A2]
    for I1, I2 in product([0,1], repeat=2):
        FE += (2*pi)**(-3)*gam2sol[I1][I2]*A[I1]*A[I2]

    for I1, I2, I3 in product([0,1], repeat=3):
        FE += (2*pi)**(-6)*(2/3/np.sqrt(3))*gam3sol.real[I1][I2][I3]*        A[I1]*A[I2]*A[I3]
        
    for I1, I2, I3, I4 in product([0,1], repeat=4):
        FE += (2*pi)**(-9)*(1/12)*        (gam41sol[I1][I2][I3][I4]+4*gam42sol[I1][I2][I3][I4]).real*        A[I1]*A[I2]*A[I3]*A[I4]
        
    return FE

def FE6sol(A1, A2, gam2sol, gam3sol, gam41sol,
           gam42sol, gam43sol, gam44sol):
    pi = 3.14159
    FE = 0
    A = [A1, A2]
    for I1, I2 in product([0,1], repeat=2):
        FE += (2*pi)**(-3)*gam2sol[I1][I2]*A[I1]*A[I2]

    for I1, I2, I3 in product([0,1], repeat=3):
        FE += (2*pi)**(-6)*(4/3/np.sqrt(6))*gam3sol.real[I1][I2][I3]*        A[I1]*A[I2]*A[I3]

    for I1, I2, I3, I4 in product([0,1], repeat=4):
        FE += (2*pi)**(-9)*(1/24)*        (gam41sol[I1][I2][I3][I4]+8*gam42sol[I1][I2][I3][I4]         +2*gam43sol[I1][I2][I3][I4]+4*gam44sol[I1][I2][I3][I4]).real*        A[I1]*A[I2]*A[I3]*A[I4]

    return FE


# In[19]:


def FE1solmin(N, FA, PHIP, Ks, CHIAB, gam41sol):
    pset=propagator.prop_set(nlam=1)
    CHIAS, CHIBS = 0, 0
    Avec0 = [10,10]
    
    gam2sol = gamma2sol(pset, N, FA, PHIP, Ks, CHIAB, CHIAS, CHIBS).real
    
    return optimize.fmin(lambda Avec: FE1sol(Avec[0], Avec[1], gam2sol, gam41sol),
                         Avec0, disp=False, full_output=1)[1]

def FE3solmin(N, FA, PHIP, Ks, CHIAB, gam3sol, gam41sol, gam42sol):
    pset=propagator.prop_set(nlam=1)
    CHIAS, CHIBS = 0, 0
    Avec0 = [10,10]
    
    gam2sol = gamma2sol(pset, N, FA, PHIP, Ks, CHIAB, CHIAS, CHIBS).real
    
    return optimize.fmin(lambda Avec: FE3sol(Avec[0], Avec[1], gam2sol,
                                      gam3sol, gam41sol, gam42sol),
                         Avec0, disp=False, full_output=1)[1]

def FE6solmin(N, FA, PHIP, Ks, CHIAB, gam3sol, gam41sol, gam42sol, gam43sol, gam44sol):
    pset=propagator.prop_set(nlam=1)
    CHIAS, CHIBS = 0, 0
    Avec0 = [10,10]
    
    gam2sol = gamma2sol(pset, N, FA, PHIP, Ks, CHIAB, CHIAS, CHIBS).real
    
    return optimize.fmin(lambda Avec: FE6sol(Avec[0], Avec[1], gam2sol,
                                      gam3sol, gam41sol, gam42sol, gam43sol, gam44sol),
                         Avec0, disp=False, full_output=1)[1]


# In[20]:


def chiootsol(N, FA, PHIP):
    pset=propagator.prop_set(nlam=5)
    wigset = wd.wigner_d_vals()

    # find critical wavemode
    CHIABs, Ks, eigvmin = gamma2sol_spinodal(pset, N, FA, PHIP)
    print('      Ks = ' + str(Ks) + ' // eigvmin = '+str(eigvmin))

    gam3sol = gamma3sol(pset, N, FA, PHIP, set3Ks(Ks)).real
    gam41sol = gamma4sol(pset, wigset, N, FA, PHIP, set4Ks1(Ks)).real
    gam42sol = gamma4sol(pset, wigset, N, FA, PHIP, set4Ks2(Ks)).real
    gam43sol = gamma4sol(pset, wigset, N, FA, PHIP, set4Ks3(Ks)).real
    gam44sol = gamma4sol(pset, wigset, N, FA, PHIP, set4Ks4(Ks)).real
        
    # calculate oot temperatures
    Avec0 = [10,10]
    FE13solsq = lambda CHIAB: (FE1solmin(N, FA, PHIP, Ks, CHIABs+CHIAB, gam41sol) -                               FE3solmin(N, FA, PHIP, Ks, CHIABs+CHIAB, gam3sol,
                                         gam41sol, gam42sol))**2
    FE36solsq = lambda CHIAB: (FE3solmin(N, FA, PHIP, Ks, CHIABs+CHIAB, gam3sol,
                                         gam41sol, gam42sol) -\
                               FE6solmin(N, FA, PHIP, Ks, CHIABs+CHIAB, gam3sol,
                                         gam41sol, gam42sol, gam43sol, gam44sol))**2

    CHI13sol = optimize.fminbound(FE13solsq, 0/N/PHIP, 100/N/PHIP)
    CHI36sol = optimize.fminbound(FE36solsq, 0/N/PHIP, 100/N/PHIP)
    
    return CHIABs, CHI13sol, CHI36sol

