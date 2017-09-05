
# coding: utf-8

# ## Vertex functions
# This code calculates the vertex functions from random-phase-approximation of copolymer melts.

# In[ ]:


from CORRcalc import s2wlc, s2inverse, s3wlc, s4wlc, norm
from itertools import product
import propagator 
import numpy as np
from scipy import optimize
import wignerD as wd


# In[ ]:


get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
# import mpld3
# mpld3.enable_notebook()


# ## Finding spinodal (minimum of $\Gamma_{2}$)

# In[4]:


def r2(N):
    return N - 0.5*(1-np.exp(-2*N))


# In[5]:


def spinodal(pset, N, FA):
    CHI = 0
    K0 = 1/np.sqrt(r2(N))
    
    KS = optimize.fmin(lambda K: np.real(gamma2(pset, N, FA, K, CHI)), K0,                      disp=False)
    
    return KS


# ## Quadratic vertex
# \begin{eqnarray}
# \Gamma_{2}(\vec{q}) \! &=& \! 
# \frac{1}{2} \left[
# -2 \chi + 
# S_{AA}^{(2)^{-1}}(\vec{q})-
# 2S_{AB}^{(2)^{-1}}(\vec{q})+
# S_{BB}^{(2)^{-1}}(\vec{q})
# \right]
# \end{eqnarray}

# In[6]:


def gamma2(pset, N, FA, K, CHI):
    s2inv = s2inverse(pset, N, FA, K)

    D = [1,-1]    # sign indicator
    G = 0
    for I0, I1 in product([0,1], repeat=2):
        G += s2inv[I0, I1]*D[I0]*D[I1]
        
    return -2*CHI + N*G


# ## Cubic Vertex
# \begin{eqnarray}
# \Gamma_{3}(\vec{q}_1,\vec{q}_2,\vec{q}_3) \! &=& \!
# -\frac{1}{3!}
# \sum_{\alpha_{1} \alpha_{2} \alpha_{3}} \! \! \!
# S_{\alpha_1 \alpha_2 \alpha_3}^{(3)}(\vec{q}_1,\vec{q}_2,\vec{q}_3) \times
# \left[ S_{\alpha_1 A}^{(2)^{-1}}(\vec{q}_1)-S_{\alpha_1 B}^{(2)^{-1}}(\vec{q}_1)  \right] \times \nonumber \\
# & &
# \hspace{.1in}
# \left[ S_{\alpha_2 A}^{(2)^{-1}}(\vec{q}_2)-S_{\alpha_2 B}^{(2)^{-1}}(\vec{q}_2)  \right]
# \left[ S_{\alpha_2 A}^{(2)^{-1}}(\vec{q}_3)-S_{\alpha_2 B}^{(2)^{-1}}(\vec{q}_3)  \right]
# \end{eqnarray}

# In[12]:


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
# 
# \begin{eqnarray}
# \Gamma_{4}(\vec{q}_1,\vec{q}_2,\vec{q}_3,\vec{q}_4) \! &=& \!
# \frac{1}{4!}
# \sum_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{3}} \! \! \!
# \gamma_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{3}}(\vec{q}_1,\vec{q}_2,\vec{q}_3,\vec{q}_4)
# \left[ S_{\alpha_1 A}^{(2)^{-1}}(\vec{q}_1)-S_{\alpha_1 B}^{(2)^{-1}}(\vec{q}_1)  \right] \times \nonumber \\
# & &
# \hspace{-1in}
# \left[ S_{\alpha_2 A}^{(2)^{-1}}(\vec{q}_2)-S_{\alpha_2 B}^{(2)^{-1}}(\vec{q}_2)  \right]
# \left[ S_{\alpha_3 A}^{(2)^{-1}}(\vec{q}_3)-S_{\alpha_3 B}^{(2)^{-1}}(\vec{q}_3)  \right]
# \left[ S_{\alpha_4 A}^{(2)^{-1}}(\vec{q}_4)-S_{\alpha_4 B}^{(2)^{-1}}(\vec{q}_4)  \right],
# \end{eqnarray}
# 
# where
# \begin{eqnarray}
# \gamma_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{3}}(\vec{q}_1,\vec{q}_2,\vec{q}_3,\vec{q}_4)
# &=&
# \sum_{\beta \gamma} \left[
# S^{(3)}_{\alpha_{1} \alpha_{2} \alpha} (\vec{q}_1, \vec{q}_2, \vec{q}_3)
# S_{\beta \gamma}^{(2)^{-1}}(\vec{q}_1+\vec{q}_2)
# S^{(3)}_{\alpha_{3} \alpha_{4} \gamma}(\vec{q}_2, 
# \vec{q}_3, \vec{q}_4) \right. \nonumber \\
# & &
# \hspace{-1in}
# + S^{(3)}_{\alpha_{1} \alpha_{2} \alpha} (\vec{q}_1, \vec{q}_3, \vec{q}_4)
# S_{\beta \gamma}^{(2)^{-1}}(\vec{q}_1+\vec{q}_3)
# S^{(3)}_{\alpha_{3} \alpha_{4} \gamma}(\vec{q}_3, 
# \vec{q}_4, \vec{q}_2) \nonumber \\
# & &
# \hspace{-1in}
# + S^{(3)}_{\alpha_{1} \alpha_{2} \alpha} (\vec{q}_2, \vec{q}_3, \vec{q}_4)
# S_{\beta \gamma}^{(2)^{-1}}(\vec{q}_2+\vec{q}_3)
# S^{(3)}_{\alpha_{3} \alpha_{4} \gamma}(\vec{q}_3, 
# \left. \vec{q}_4, \vec{q}_1) \right]
# - S^{(4)}_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{3}}(\vec{q}_1,\vec{q}_2,\vec{q}_3,\vec{q}_4)
# \end{eqnarray}

# In[13]:


def gamma4(pset, wigset, N, FA, Ks):
    K1, K2, K3, K4 = Ks
    if not (abs(norm(K1)-norm(K2)) < 1e-5
            and abs(norm(K2)-norm(K3)) < 1e-5
            and abs(norm(K3)-norm(K4)) < 1e-5):
        print(K1, K2, K3, K4)
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
            G4[a1][a2][a3][a4] +=                 s31[a1][a2][I0]*s31[a3][a4][I1]*s21inv[I0][I1] +                 s32[a1][a4][I0]*s32[a2][a3][I1]*s22inv[I0][I1] +                 s33[a1][a3][I0]*s33[a2][a4][I1]*s23inv[I0][I1]
    G4 -= s4
    
    val = 0
    for I0, I1, I2, I3 in product([0,1], repeat=4):
        val += G4[I0][I1][I2][I3] *                (s2inv[I0][0] - s2inv[I0][1])*                (s2inv[I1][0] - s2inv[I1][1])*                (s2inv[I2][0] - s2inv[I2][1])*                (s2inv[I3][0] - s2inv[I3][1])
                
    return val*(N**3)


# ## Quadratic vertex for copolymer solutions

# In[100]:


def gamma2sol(pset, N, FA, PHIP, K, CHIAB, CHIAS, CHIBS):    
    CHIMAT = np.zeros((2, 2))
    CHIMAT[0][0] = CHIAS
    CHIMAT[0][1] = 0.5*(-CHIAB+CHIAS+CHIBS)
    CHIMAT[1][0] = CHIMAT[0][1]
    CHIMAT[1][1] = CHIBS
    s2inv = s2inverse(pset, N, FA, K)
    
    G = np.zeros((2, 2), dtype=type(1+1j))
    for I0, I1 in product([0,1], repeat=2):
        G[I0][I1] = -2*CHIMAT[I0][I1] + N*s2inv[I0][I1]/PHIP + 1/(1-PHIP)

    # where A is the inverse matrix of Fredrickson-Leibler's 1989 solution paper
    A = np.array([[PHIP, FA], [-PHIP, 1-FA]], dtype=type(1+1j))
    return np.matmul(np.matmul(A.T, G), A)


# ## Cubic vertex for copolymer solutions

# In[132]:


def Qsol(pset, N, FA, PHIP, K):
    s2inv = s2inverse(pset, N, FA, norm(K))
    
    Q = np.zeros((2, 2), dtype = type(1+1j))
    Q[0][0] = s2inv[0][0] - s2inv[0][1]
    Q[0][1] = s2inv[0][1] - s2inv[1][1]
    Q[1][0] = (1/PHIP)*(FA*s2inv[0][0] + (1-FA)*s2inv[0][1])
    Q[1][1] = (1/PHIP)*(FA*s2inv[0][1] + (1-FA)*s2inv[1][1])

    return Q
    
def gamma3sol(pset, N, FA, PHIP, Ks):
    K1, K2, K3 = Ks
    if norm(K1+K2+K3) >= 1e-10:
        raise('Qs must add up to zero')
        
#     if not (abs(norm(K1)-norm(K2)) < 1e-5 \
#         and abs(norm(K2)-norm(K3)) < 1e-5):
#         raise('Qs must have same length')
    
    s3 = s3wlc(pset, N, FA, Ks)
    Q1 = Qsol(pset, N, FA, PHIP, K1)
    Q2 = Qsol(pset, N, FA, PHIP, K2)
    Q3 = Qsol(pset, N, FA, PHIP, K3)
    
    G3 = np.zeros((2, 2, 2), dtype=type(1+1j))
    for I0, I1, I2 in product([0,1], repeat=3):
        for J0, J1, J2 in product([0,1], repeat=3):
            G3[I0][I1][I2] -= PHIP* s3[J0][J1][J2]*                                    Q1[I0][J0]* Q2[I1][J1]* Q3[I2][J2]

    return G3*(N**2)


# ## Quartic vertex for copolymer solutions

# In[143]:


def gamma4sol(pset, wigset, N, FA, PHIP, Ks):
    K1, K2, K3, K4 = Ks
#     if not (norm(K1) == norm(K2) == norm(K3) == norm(K4)):
#         raise('Qs must have same length')
    
    K = norm(K1)
    K12 = norm(K1+K2)
    K13 = norm(K1+K3)
    K14 = norm(K1+K4)
    
    s4 = s4wlc(pset, wigset, N, FA, Ks)
    s31 = s3wlc(pset, N, FA, [K1, K2, -K1-K2])
    s32 = s3wlc(pset, N, FA, [K1, K3, -K1-K3])
    s33 = s3wlc(pset, N, FA, [K1, K4, -K1-K4])

    s21inv = s2inverse(pset, N, FA, K12)
    s22inv = s2inverse(pset, N, FA, K13)
    s23inv = s2inverse(pset, N, FA, K14)
    
    Q1 = Qsol(pset, N, FA, PHIP, K1)
    Q2 = Qsol(pset, N, FA, PHIP, K2)
    Q3 = Qsol(pset, N, FA, PHIP, K3)
    Q4 = Qsol(pset, N, FA, PHIP, K4)

    G = np.zeros((2,2,2,2),dtype=type(1+1j))
    for a1, a2, a3, a4 in product([0,1], repeat=4):
        for I0, I1 in product([0,1], repeat=2):
            G[a1][a2][a3][a4] +=             s31[a1][a2][I0]*s31[a3][a4][I1]*s21inv[I0][I1] +             s32[a1][a4][I0]*s32[a2][a3][I1]*s22inv[I0][I1] +             s33[a1][a3][I0]*s33[a2][a4][I1]*s23inv[I0][I1]
    G -= s4
    
    G4 = np.zeros((2, 2, 2, 2), dtype=type(1+1j))
    for I0, I1, I2, I3 in product([0,1], repeat=4):
        for J0, J1, J2, J3 in product([0,1], repeat=4):
            G4[I0][I1][I2][I3] = PHIP* G[J0][J1][J2][J3]*                            Q1[I0][J0]*Q2[I1][J1]*Q3[I2][J2]*Q4[I3][J3]

    return G4*(N**3)

