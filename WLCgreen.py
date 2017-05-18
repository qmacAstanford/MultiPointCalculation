
# coding: utf-8

# # Eigenvalues and Summand for WLC

# This is a python implementation of eigenvalue and summand evaluation for the worm like chain propagator based on the 2007 paper "End-to-End distribution for a wormlike chain in arbitrary dimensions" by Shafigh Mehraeen, Bariz Sudhanshu, Elena Koslover, and Andrew Spakowitz.
# This code was translated into python from matlab by Quinn in 2017.

# In[1]:

import numpy as np
from numba import jit


# $
# a_j^\mu=\sqrt{\frac{(j-\mu)(j+\mu+D-3)}{(2j+D-2)(2j+D-4)}}
# $

# In[2]:

@jit(nopython=True)
def get_a(lam,mu,d=3):
    return np.sqrt((lam-mu)*(lam+mu+d-3)/                   float((2*lam+d-2)*(2*lam+d-4)))


# ## Eigenvalues

# In[2]:

# eigvals = wlc.wlc_eigvals(k,ORDEig,mu,d=3)


# ### Intermediate K

# We would like to calculate the eigenvalues denoted $\epsilon_l$ in the paper.  These values depend on the z-component quantum number $\mu$ and the magnitude of the Fourier frequency $K$.

# In the intermediate K regime we use the method in appendix A.2. We find the eigenvalues of the tri-diagonal matrix $J$ where the diagonal elements are 
# $$
# J_{j,j}=(\mu+j-1)(\mu+j+D-3)
# $$
# And off diagonal elements by
# $$
# J_{j+1,j}=J_{j,j+1}=-ia_{\mu+j}^{\mu}K
# $$
# only the eigenvalues ORDEig eigenvalues with the largest real part are returned.

# In[4]:

@jit(nopython=True)
def make_E(n,mu,K,d=3):
    E=np.zeros((n,n),dtype=type(1+1j))
    left=(0+0j)*np.NaN
    for row in range(0,n):
        j=row+1
        if K<=1:
            right=-1j*K*get_a(j+mu,mu,d) # for above diagonal
            if row>0:
                left=-1j*K*get_a(j+mu-1,mu,d) # for below diagonal
            diag=(mu+j-1)*(mu+j+d-3) # diagonal element
            if row==0:
                E[row,0:2]=[diag,right]
            elif row==n-1:
                E[row,row-1:row+1]=[left,diag]
            else:
                E[row,row-1:row+2]=[left,diag,right]
        else:
            right=-1j*get_a(j+mu,mu,d) # for above diagonal
            if row>0:
                left=-1j*get_a(j+mu-1,mu,d) # for below diagonal
            diag=(mu+j-1)*(mu+j+d-3.0)/K # diagonal element
            if row==0:
                E[row,0:2]=[diag,right]
            elif row==n-1:
                E[row,row-1:row+1]=[left,diag]
            else:
                E[row,row-1:row+2]=[left,diag,right]
    return E


# In[5]:

# Returns top ORDEig eigenvalus using Intermediate K methode
#@jit #(nopython=True)
def IntermediateKEigenValues(K,ORDEig,mu,d=3):
    N=int(2*(np.ceil(ORDEig/2.0))) # Eigenvalues come in pairs
    
    n=int(4*N) # Calculate extra 
    E=make_E(n,mu,K,d)
    #import scipy.sparse as sp
    #from scipy.sparse import linalg as sp_linalg
    #dmtrx = sp.dia_matrix.todia(E)
    #evals = sp_linalg.eigs(dmtrx,N,which='LR',return_eigenvectors=False)
    evals = -1*np.linalg.eigvals(E)    
    evals = np.sort(evals)[::-1]
    
    # Fix sorting
    for ii in range(0,N,2):
        evals[ii]=np.real(evals[ii])+1j*abs(np.imag(evals[ii]))
        evals[ii+1]=np.real(evals[ii+1])-1j*abs(np.imag(evals[ii+1]))
    
    if K>1:
        evals=evals*K
    
    out=np.zeros(ORDEig,'complex')*np.NaN
    out[mu:]=evals[0:ORDEig-mu]
    return out
    


# ### Large K

# In[6]:

# implementation equivalent to table II fore E
def Epsilon(l,d,alpha,mu):
    beta=-np.sqrt(2)/4*(1+1j);
    m=mu+(d-3)/2.0;  # Q.J.M. changes this line 8/1/15
    n=2*l+m+1;  # also know as s

    epsilon_0=(-1.0/2/beta)**(-1)*(n/2.0);
    epsilon_1=(-1.0/2/beta)**( 0)*(-1.0/8*(n**2+3-3*m**2)-m*(m+1));
    epsilon_2=(-1.0/2/beta)**( 1)*(-1.0/2**5*n*(n**2+3-9*m**2));
    epsilon_3=(-1.0/2/beta)**( 2)*(-1.0/2**8*(5*n**4+34*n**2+9)-(102*n**2+42)*m**2+33*m**4);
    epsilon_4=(-1.0/2/beta)**( 3)*(-1.0/2**11*n*(33*n**4+410*n**2+405)-(1230*n**2+1722)*m**2                                        +813*m**4);
    epsilon_5=(-1.0/2/beta)**( 4)*(-1.0/2**12*9*(7*n**6+140*n**4+327*n**2                                        +54-(420*n**4+1350*n**2+286)*m**2                                        +(495*n**2+314)*m**4-82*m**6));

    value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha**2+           epsilon_4*alpha**3+epsilon_5*alpha**4;
    return value


# Large K eigenvalues are implemented via
# $
#    \epsilon_{2l}=iK-\mu(\mu+D-2)-\frac{1}{\alpha}\sum_{n=0}^{\infty}E_n^{(l)}\alpha^n
# $

# In[7]:

def LargeKEigenValues(k,ORDEig,mu,d=3):
    Eig=np.zeros(ORDEig+1,dtype='complex')*np.NaN
    for l in range(0,ORDEig,2):
        alpha=1.0/np.sqrt(8*k)
        Eig[l]=1j*k-mu*(mu+d-2)-Epsilon(l/2.0,d,alpha,mu)
        Eig[l+1]=np.conj(Eig[l])
        
    out=np.zeros(ORDEig,'complex')*np.NaN
    out[mu:]=Eig[0:ORDEig-mu]
    return out
    


# In[8]:

# Cutoff based on graph
def wlc_eigvals(k,ORDEig,mu,d=3):
    cutoff=800
    if abs(k)<cutoff:
        Eig=IntermediateKEigenValues(k,ORDEig,mu,d)
    else:
        Eig=LargeKEigenValues(k,ORDEig,mu,d)
    return Eig


# ## GLM

# ### G

# In[3]:

# sequence of integers from first to last
# for example seq(2,5)=2,3,4,5
@jit(nopython=True)
def seq(first,last):
    if first<=last:
        return np.arange(first,last+1)
    else:
        return np.arange(first,last-1,-1)


# $$
# P_{\lambda}=p+\lambda\left(\lambda+d-2\right)
# $$

# In[4]:

@jit(nopython=True)
def get_P(lam,p,d=3):
    return p+lam*(lam+d-2)


# $$
# j_{\lambda}^{\mu\left(+\right)}=P_{\lambda}+\left(a_{\lambda+1}^{\mu}K\right)^{2}\frac{1}{j_{\lambda+1}^{\mu\left(+\right)}}
# $$

# In[15]:

# Calculate j(+), output indexed by lambda size: lamMax+1
@jit(nopython=True)
def get_jp(p,mu,K,lamMax,d=3):
    j=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if K<=1.0:
        j[lamMax]=get_P(lamMax,p,d=d) # Just guessing
        for lam in seq(lamMax-1,mu):
            j[lam] = get_P(lam,p,d=d)+(get_a(lam+1,mu,d=d)*K)**2 /j[lam+1]
        return j
    else:
        j[lamMax]=get_P(lamMax,p,d=d)/K # Just guessing
        for lam in seq(lamMax-1,mu):
            j[lam] = get_P(lam,p,d=d)/K+(get_a(lam+1,mu,d=d))**2 /j[lam+1]
        j=j*K
        return j


# $$
# j_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# P_{\lambda}+\left(a_{\lambda}^{\mu}K\right)^{2}\frac{1}{j_{\lambda-1}^{\mu\left(-\right)}} & \lambda\geq\mu+1\\
# P_{\lambda} & \lambda=\mu
# \end{cases}
# $$

# In[6]:

@jit(nopython=True)
def get_jm(p,mu,K,lamMax,d=3):
    jm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if K<=1.0:
        jm[mu]=get_P(mu,p,d=d) 
        for lam in seq(mu+1,lamMax):
            jm[lam]=get_P(lam,p,d=d) + (get_a(lam,mu,d=d)*K)**2 /jm[lam-1]
        return jm
    else:
        jm[mu]=get_P(mu,p,d=d)/K 
        for lam in seq(mu+1,lamMax):
            jm[lam]=get_P(lam,p,d=d)/K + (get_a(lam,mu,d=d))**2 /jm[lam-1]
        jm=jm*K
        return jm


# $$
# \frac{\partial}{\partial p}j_{\lambda}^{\mu\left(+\right)}=1-\left(a_{\lambda+1}^{\mu}K\right)^{2}\frac{1}{\left(j_{\lambda+1}^{\mu\left(+\right)}\right)^{2}}\frac{\partial}{\partial p}j_{\lambda+1}^{\mu\left(+\right)}
# $$

# In[7]:

@jit(nopython=True)
def get_djp(p,mu,K,lamMax,jp,d=3):
    djp=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if K<=1.0:
        djp[lamMax]= 1.0 # Just guessing
        for lam in seq(lamMax-1,mu):
            djp[lam] = 1.0 - (get_a(lam+1,mu,d=d)*K)**2 * djp[lam+1]/(jp[lam+1]**2)
        return djp 
    else:
        djp[lamMax]= 1.0/K # Just guessing
        for lam in seq(lamMax-1,mu):
            djp[lam] = 1.0/K - (get_a(lam+1,mu,d=d)*K)**2 * djp[lam+1]/(jp[lam+1]**2)
        djp=djp*K
        return djp 


# $$
# \frac{\partial}{\partial p}j_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# 1-\left(a_{\lambda}^{\mu}K\right)^{2}\frac{1}{\left(j_{\lambda-1}^{\mu\left(-\right)}\right)^{2}}\frac{\partial}{\partial p}j_{\lambda-1}^{\mu\left(-\right)} & \lambda\geq\mu+1\\
# 1 & \lambda=\mu
# \end{cases}
# $$

# In[8]:

@jit(nopython=True)
def get_djm(p,mu,K,lamMax,jm,d=3):
    djm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if K<=1.0:
        djm[mu]= 1.0
        for lam in seq(mu+1,lamMax):
            djm[lam] = 1.0 - (get_a(lam,mu,d=d)*K)**2 * djm[lam-1]/(jm[lam-1]**2)
        return djm
    else:
        djm[mu]= 1.0/K
        for lam in seq(mu+1,lamMax):
            djm[lam] = 1.0/K - (get_a(lam,mu,d=d)*K)**2 * djm[lam-1]/(jm[lam-1]**2)
        djm=djm*K
        return djm


# $$
# \partial_{p}w_{\lambda}^{\mu\left(\pm\right)}=\frac{-1}{\left(j_{\lambda}^{\mu\left(\pm\right)}\right)^{2}}\partial_{p}j_{\lambda}^{\mu\left(\pm\right)}
# $$

# In[ ]:

@jit(nopython=True)
def get_dw(j,dj,mu):
    dw=-dj/(j**2)
    #dw=j*np.NaN
    #dw[mu:]=-dj[mu:]/(j[mu:]**2)
    return dw


# $$
# W_{\lambda}^{\mu}=\frac{1}{\left(a_{\lambda}^{\mu}K\right)^{2}w_{\lambda-1}^{\mu\left(-\right)}+P+\left(a_{\lambda+1}^{\mu}K\right)^{2}w_{\lambda+1}^{\mu\left(+\right)}}
# $$

# In[82]:

@jit(nopython=True)
def get_W(p,mu,K,lamMax,jp,jm,d=3):
    W=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    lam=mu
    W[lam]=get_P(lam,p,d=d)+(get_a(lam+1,mu,d=d)*K)**2 / jp[lam+1]
    for lam in seq(mu+1,lamMax-1):
        W[lam]=(get_a(lam,mu,d=d)*K)**2 /jm[lam-1] +              get_P(lam,p,d=d)+              (get_a(lam+1,mu,d=d)*K)**2 / jp[lam+1]
    return (1/W)


# $$
# \partial_{p}W_{\lambda}^{\mu}=-\left(W_{\lambda}^{\mu}\right)^{2}\left(\left(a_{\lambda}^{\mu}K\right)^{2}\partial_{p}w_{\lambda-1}^{\mu\left(-\right)}+1+\left(a_{\lambda+1}^{\mu}K\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\right)}\right)
# $$

# In[83]:

@jit(nopython=True)
def get_dW(p,mu,K,lamMax,dwp,dwm,W,d=3):
    dW=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    lam=mu
    dW[lam]=1.0 + (get_a(lam+1,mu,d=d)*K)**2 * dwp[lam+1]
    for lam in seq(mu+1,lamMax-1):
        dW[lam]=(get_a(lam,mu,d=d)*K)**2 *dwm[lam-1] +              1.0 +              (get_a(lam+1,mu,d=d)*K)**2 * dwp[lam+1]
    return (-dW*(W**2))


# \[
# \mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\begin{cases}
# W_{\lambda_{0}}^{\mu} & \lambda_{0}=\lambda\\
# W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\right)} & \lambda_{0}<\lambda\\
# W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)} & \lambda_{0}>\lambda
# \end{cases}
# \]

# In[84]:

@jit(nopython=True)
def get_G(p,lam,lam0,mu,K,jm,jp,W,d=3):
    if lam==lam0:
        return W[lam0]
    elif lam0<lam:
        out=W[lam0]
        for lamt in seq(lam0+1,lam):
            out=out*1j*K*get_a(lamt,mu,d=d)/jp[lamt]
        return out
    else:
        out=W[lam0]
        for lamt in seq(lam+1,lam0):
            out=out*1j*K*get_a(lamt,mu,d=d)/jm[lamt-1]
        return out


# $$
# \partial_{p}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\begin{cases}
# \partial_{p}W_{\lambda_{0}}^{\mu} & \lambda_{0}=\lambda\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}\frac{\partial_{p}w_{\tilde{\lambda}}^{\mu\left(+\right)}}{w_{\tilde{\lambda}}^{\mu\left(+\right)}}\right)\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\right)} & \lambda_{0}<\lambda\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}\right)\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)} & \lambda_{0}>\lambda
# \end{cases}
# $$

# In[85]:

@jit(nopython=True)
def get_dG(p,lam,lam0,mu,K,jm,jp,W,dW,dwm,dwp,d=3):
    if lam==lam0:
        return dW[lam0]
    elif lam0<lam:
        Sum=0.0
        Prod=1.0
        for lamt in seq(lam0+1,lam):
            Prod=Prod*1j*K*get_a(lamt,mu,d=d)/jp[lamt]
            Sum=Sum+dwp[lamt] * jp[lamt]
        return ((dW[lam0] + W[lam0]*Sum)*Prod)
    else:
        Sum=0.0
        Prod=1.0
        for lamt in seq(lam+1,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d=d)/jm[lamt-1]
            Sum=Sum+dwm[lamt-1] * jm[lamt-1]
        
        return ((dW[lam0] + W[lam0]*Sum)*Prod)     


# ### Lim p -> pole

# \[
# X_{1}=\left(a_{\lambda_{0}}^{\mu}K\right)^{2}w_{\lambda_{0}-1}^{\mu\left(-\right)}+P_{\lambda_{0}}+\left(a_{\lambda_{0}+1}^{\mu}K\right)^{2}w_{\lambda_{0}+1}^{\mu\left(+\right)}
# \]

# In[20]:

@jit(nopython=True)
def X1(p,lam0,mu,K,jp,jm,d=3):
    if lam0 == mu:
        out=  0             + get_P(lam0,p,d=d)             +(get_a(lam0+1,mu,d=d)*K)**2  * 1/jp[lam0+1]
    else:
        out=  (get_a(lam0,mu,d=d)*K)**2 * 1/jm[lam0-1]            + get_P(lam0,p,d=d)             +(get_a(lam0+1,mu,d=d)*K)**2  * 1/jp[lam0+1]
    return out


# \[
# \frac{\partial}{\partial p}X_{1}=\left(a_{\lambda_{0}}^{\mu}K\right)^{2}\frac{\partial}{\partial p}w_{\lambda_{0}-1}^{\mu\left(-\right)}+1+\left(a_{\lambda_{0}+1}^{\mu}K\right)^{2}\frac{\partial}{\partial p}w_{\lambda_{0}+1}^{\mu\left(+\right)}
# \]

# In[21]:

@jit(nopython=True)
def X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3):
    if lam0 == mu:
        X1=  0             + 1             +(get_a(lam0+1,mu,d=d)*K)**2  * get_dw(jp[lam0+1],djp[lam0+1],mu)           
    else:
        X1=  (get_a(lam0,mu,d=d)*K)**2 * get_dw(jm[lam0-1],djm[lam0-1],mu)            + 1             +(get_a(lam0+1,mu,d=d)*K)**2  * get_dw(jp[lam0+1],djp[lam0+1],mu)    
    return X1


# \[
# X^{\left(+\right)}=\prod_{n=1}^{\lambda-\lambda_{0}}j_{\lambda_{0}+n}^{\mu\left(+\right)}
# \]

# In[22]:

@jit(nopython=True)
def Xp(lam,lam0,jp):
    prod=1
    for n in seq(1,lam-lam0):
        prod=prod*jp[lam0+n]
    return prod


# \[
# X^{\left(-\right)}=\prod_{n=1}^{\lambda_{0}-\lambda}j_{\lambda_{0}-n}^{\mu\left(-\right)}
# \]

# In[23]:

@jit(nopython=True)
def Xm(lam,lam0,jm):
    prod=1
    for n in seq(1,lam0-lam):
        prod=prod*jm[lam0-n]
    return prod


# \[
# c^{\left(+\right)}=\prod_{n=1}^{\lambda-\lambda_{0}}iKa_{\lambda_{0}+n}^{\mu}
# \]

# In[24]:

@jit(nopython=True)
def get_cp(lam,lam0,mu,K,d=3):
    prod=1
    for n in seq(1,lam-lam0):
        prod=prod*1j*K*get_a(lam0+n,mu,d)
    return prod


# \[
# c^{\left(-\right)}=\prod_{n=1}^{\lambda_{0}-\lambda}iKa_{\lambda_{0}+1-n}^{\mu}
# \]

# In[25]:

@jit(nopython=True)
def get_cm(lam,lam0,mu,K,d=3):
    prod=1
    for n in seq(1,lam0-lam):
        prod=prod*1j*K*get_a(lam0+1-n,mu,d)
    return prod


# \[
# \frac{\partial}{\partial p}X^{\left(+\right)}=X^{\left(+\right)}\left(\sum_{n=1}^{\lambda-\lambda_{0}}\frac{\frac{\partial}{\partial p}j_{\lambda_{0}+n}^{\mu\left(+\right)}}{j_{\lambda_{0}+n}^{\mu\left(+\right)}}\right)
# \]

# In[26]:

@jit(nopython=True)
def Xpprime(lam,lam0,jp,djp,Xp):
    out=0.0
    for n in seq(1,lam-lam0):
        out=out+(djp[lam0+n]/jp[lam0+n])
    out=out*Xp
    return out


# \[
# \frac{\partial}{\partial p}X^{\left(-\right)}=X^{\left(-\right)}\left(\sum_{n=1}^{\lambda_{0}-\lambda}\frac{\frac{\partial}{\partial p}j_{\lambda_{0}-n}^{\mu\left(-\right)}}{j_{\lambda_{0}-n}^{\mu\left(-\right)}}\right)
# \]

# In[27]:

@jit(nopython=True)
def Xmprime(lam,lam0,jm,djm,Xm):
    out=0.0
    for n in seq(1,lam0-lam):
        out=out+(djm[lam0-n]/jm[lam0-n])
    out=out*Xm
    return out


# \[
# \frac{1}{\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(p\right)}=\begin{cases}
# X_{1} & \lambda_{0}=\lambda\\
# \frac{X_{1}X^{\left(+\right)}}{c^{\left(+\right)}} & \lambda_{0}<\lambda\\
# \frac{X_{1}X^{\left(-\right)}}{c^{\left(-\right)}} & \lambda_{0}>\lambda
# \end{cases}
# \]

# In[28]:

@jit(nopython=True)
def InvG(p,lam,lam0,mu,K,jp,jm,d=3):
    if min(lam,lam0)<mu:
        return np.NaN
    if lam0==lam:
        return X1(p,lam0,mu,K,jp,jm)
    if lam0<lam:
        return (X1(p,lam0,mu,K,jp,jm)*Xp(lam,lam0,jp)/get_cp(lam,lam0,mu,K,d))
    if lam0>lam:
        return (X1(p,lam0,mu,K,jp,jm)*Xm(lam,lam0,jm)/get_cm(lam,lam0,mu,K,d))


# \[
# \frac{\partial}{\partial p}\left(\frac{1}{\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(p\right)}\right)=\begin{cases}
# \partial_{p}X_{1} & \lambda_{0}=\lambda\\
# \frac{X_{1}\partial_{p}X^{\left(+\right)}+X^{\left(+\right)}\partial_{p}X_{1}}{c^{\left(+\right)}} & \lambda_{0}<\lambda\\
# \frac{X_{1}\partial_{p}X^{\left(-\right)}+X^{\left(-\right)}\partial_{p}X_{1}}{c^{\left(-\right)}} & \lambda_{0}>\lambda
# \end{cases}
# \]

# In[29]:

@jit(nopython=True)
def InvGprimeGeneral(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3):
    if min(lam,lam0)<mu:
        return np.NaN
    if lam0==lam:
        return X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)
    if lam0<lam:
        xp=Xp(lam,lam0,jp)
        return ((X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*xp                 +X1(p,lam0,mu,K,jp,jm)*Xpprime(lam,lam0,jp,djp,xp))                /get_cp(lam,lam0,mu,K,d))
    if lam0>lam:
        xm=Xm(lam,lam0,jm)
        return ((X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*xm                 +X1(p,lam0,mu,K,jp,jm)*Xmprime(lam,lam0,jm,djm,xm))                /get_cm(lam,lam0,mu,K,d))


# \[
# \frac{\partial}{\partial p}\left(\frac{1}{\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(p\right)}\right)=\begin{cases}
# \partial_{p}X_{1} & \lambda_{0}=\lambda\\
# \frac{X^{\left(+\right)}\partial_{p}X_{1}}{c^{\left(+\right)}} & \lambda_{0}<\lambda\\
# \frac{X^{\left(-\right)}\partial_{p}X_{1}}{c^{\left(-\right)}} & \lambda_{0}>\lambda
# \end{cases}
# \]
# because $X_1=0$ at the poles, $\partial_p X^\left{ \pm \right}$ does not generally diverge at the poles.

# In[30]:

@jit(nopython=True)
def InvGprime(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3):
    if min(lam,lam0)<mu:
        return np.NaN
    tol=0.00001
    if abs(X1(p,lam0,mu,K,jp,jm)) > tol:
        print("Error in InvGpriem. X1 should be zero at pole")
        return np.NaN
    
    if lam0==lam:
        return X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)
    if lam0<lam:
        return (X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*Xp(lam,lam0,jp)                /get_cp(lam,lam0,mu,K,d))
    if lam0>lam:
        return (X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*Xm(lam,lam0,jm)                /get_cm(lam,lam0,mu,K,d))


# ### Limit p -> 0

# $$
# \lim_{p\to0}j_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# p & \lambda=\mu=0\\
# \left(a_{\lambda}^{\mu}K\right)^{2}/p & \lambda=1,\ \mu=0\\
# \lambda\left(\lambda+1\right) & \lambda=2,\ \mu=0\\
# \lambda\left(\lambda+1\right) & \lambda=\mu>0\\
# \lambda\left(\lambda+1\right)+\left(a_{\lambda}^{\mu}K\right)^{2}\frac{1}{\lim_{p\to0}j_{\lambda-1}^{\mu\left(-\right)}} & \mathrm{otherwise}
# \end{cases}
# $$

# In[86]:

@jit(nopython=True)
def get_jm_zero(mu,K,lamMax,d=3):
    jm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu>0:    
        jm[mu]=get_P(mu,0,d) 
        for lam in seq(mu+1,lamMax):
            jm[lam]=get_P(lam,0,d) + (get_a(lam,mu,d)*K)**2 /jm[lam-1]
        return jm
    elif mu==0:
        jm[0]=0
        jm[1]=np.Inf
        jm[2]=6.0
        for lam in seq(3,lamMax):
            jm[lam]=get_P(lam,0,d) + (get_a(lam,mu,d)*K)**2 /jm[lam-1]
        return jm


# $$
# \lim_{p\to0}w_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# 1/p & \lambda=\mu=0\\
# p/\left(a_{\lambda}^{\mu}K\right)^{2} & \lambda=1,\ \mu=0\\
# 1/\lim_{p\to0}j_{\lambda}^{\mu\left(-\right)} & \mathrm{\mathrm{otherwise}}
# \end{cases}
# $$

# In[87]:

@jit(nopython=True)
def get_wm_zero(jm,lamMax):
    wm=np.zeros(lamMax+1,dtype=type(1+1j))
    wm[0]=(1+0j)*np.Inf
    wm[2:lamMax+1]=1.0/jm[2:lamMax+1]


# $$
# \lim_{p\to0}\frac{\partial}{\partial p}j_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# 1 & \lambda=\mu=0\\
# -\left(a_{\lambda}^{\mu}K\right)^{2}/p^{2} & \lambda=1,\ \mu=0\\
# 1+\frac{\left(a_{\lambda}^{\mu}\right)^{2}}{\left(a_{\lambda-1}^{\mu}\right)^{2}} & \lambda=2,\ \mu=0\\
# 1 & \lambda=\mu>0\\
# 1-\left(a_{\lambda}^{\mu}K\right)^{2}\frac{1}{j_{\lambda-1}^{\mu\left(-\right)}{}^{2}}\frac{\partial}{\partial p}j_{\lambda-1}^{\mu\left(-\right)} & \mathrm{otherwise}
# \end{cases}
# $$

# In[88]:

@jit(nopython=True)
def get_djm_zero(mu,K,lamMax,jm,d=3):
    djm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu>0:
        djm[mu]= 1.0
        for lam in seq(mu+1,lamMax):
            djm[lam] = 1.0 - (get_a(lam,mu,d)*K)**2 * djm[lam-1]/(jm[lam-1]**2)
        return djm    
    elif mu==0:
        djm[0]=1.0
        djm[1]=(1+0j)*np.Inf
        djm[2]=1.0+(get_a(2,mu,d)**2)/(get_a(1,mu,d)**2)
        for lam in seq(3,lamMax):
            djm[lam] = 1.0 - (get_a(lam,mu,d)*K)**2 * djm[lam-1]/(jm[lam-1]**2)
        return djm


# $$
# \lim_{p\to0}\partial_{p}w_{\lambda}^{\mu\left(-\right)}=\begin{cases}
# -\frac{1}{p^{2}} & \lambda=\mu=0\\
# \frac{1}{\left(a_{\lambda}^{\mu}K\right)^{2}} & \lambda=1,\ \mu=0\\
# \frac{-1}{\left(j_{\lambda}^{\mu\left(-\right)}\right)^{2}}\partial_{p}j_{\lambda}^{\mu\left(-\right)} & \mathrm{otherwise}
# \end{cases}
# $$

# In[90]:

@jit(nopython=True)
def get_dwm_zero(jm,djm,K,lamMax,mu,d=3):
    dwm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu>0:
        dwm[mu:]=-1.0*djm[mu:]/(jm[mu:]**2)
    elif mu==0:
        dwm[0]=(-1.0+0j)*np.Inf
        dwm[1]=1/((get_a(1,0,d)*K)**2)
        dwm[2:lamMax+1]=-1.0*djm[2:lamMax+1]/(jm[2:lamMax+1]**2)
    return dwm


# $$
# \lim_{p\to0}W_{\lambda}^{\mu}=\begin{cases}
# \frac{1}{\lambda\left(\lambda+1\right)+\left(a_{\lambda+1}^{\mu}K\right)^{2}w_{\lambda+1}^{\mu\left(+\right)}} & \lambda=0,2,\ \mu=0\\
# \frac{p}{\left(a_{\lambda}^{\mu}K\right)^{2}} & \lambda=1,\ \mu=0\\
# \frac{1}{\left(a_{\lambda}^{\mu}K\right)^{2}w_{\lambda-1}^{\mu\left(-\right)}+\lambda\left(\lambda+1\right)+\left(a_{\lambda+1}^{\mu}K\right)^{2}w_{\lambda+1}^{\mu\left(+\right)}} & \mathrm{otherwise}
# \end{cases}
# $$

# In[91]:

@jit(nopython=True)
def get_W_zero(mu,K,lamMax,jp,jm,d=3):
    W=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu>0:
        lam=mu
        W[lam]=get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1]
        for lam in seq(mu+1,lamMax-1):
            W[lam]=(get_a(lam,mu,d)*K)**2 /jm[lam-1] +                  get_P(lam,0,d)+                  (get_a(lam+1,mu,d)*K)**2 / jp[lam+1]
        return (1/W)
    elif mu==0:
        lam=0
        W[0]=1.0/(get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        W[1]=0.0
        lam=2
        W[2]=1.0/(get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        for lam in seq(3,lamMax-1):
            W[lam]=1.0/((get_a(lam,mu,d)*K)**2 /jm[lam-1] +                  get_P(lam,0,d)+                  (get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        return W 


# $$
# \lim_{p\to0}\partial_{p}W_{\lambda}^{\mu}=\begin{cases}
# -\left(W_{\lambda}^{\mu}\right)^{2}\left(1+\left(a_{\lambda+1}^{\mu}K\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\right)}\right) & \lambda=0,\ \mu=0\\
# \frac{1}{\left(a_{\lambda}^{\mu}K\right)^{2}} & \lambda=1,\ \mu=0\\
# -\left(W_{\lambda}^{\mu}\right)^{2}\left(\left(a_{\lambda}^{\mu}K\right)^{2}\partial_{p}w_{\lambda-1}^{\mu\left(-\right)}+1+\left(a_{\lambda+1}^{\mu}K\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\right)}\right) & \mathrm{otherwise}
# \end{cases}
# $$

# In[92]:

@jit(nopython=True)
def get_dW_zero(mu,K,lamMax,dwp,dwm,W,d=3):
    dW=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu>0:
        lam=mu
        dW[lam]=1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1]
        for lam in seq(mu+1,lamMax-1):
            dW[lam]=(get_a(lam,mu,d)*K)**2 *dwm[lam-1] +                  1.0 +                  (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1]
        return (-dW*(W**2))
    elif mu==0:
        lam=0
        dW[0]=-(W[lam]**2)*(1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1])
        dW[1]=(get_a(1,0,d)*K)**(-2)
        for lam in seq(2,lamMax-1):
            dW[lam]=-(W[lam]**2)*((get_a(lam,mu,d)*K)**2 *dwm[lam-1] +                                   1.0 +                                   (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1])
        return dW


# $$
# \lim_{p\to0}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\begin{cases}
# W_{\lambda_{0}}^{\mu} & \lambda_{0}=\lambda\\
# W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\right)} & \lambda_{0}<\lambda\\
# \frac{i}{a_{1}^{0}K} & \lambda_{0}=1,\ \lambda=\mu=0\\
# W_{\lambda_{0}}^{\mu}\left(-\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\right) & \lambda_{0}=2,\ \lambda=\mu=0\\
# W_{\lambda_{0}}^{\mu}\left(-\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\right)\left(\prod_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}\right) & \lambda_{0}>2,\ \lambda=\mu=0\\
# W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)} & \lambda_{0}>\lambda>0
# \end{cases}
# $$

# In[37]:

@jit(nopython=True)
def get_G_zero(lam,lam0,mu,K,jm,jp,W,d=3):
    if lam==lam0:
        return W[lam0]
    elif lam0<lam:
        out=W[lam0]
        for lamt in seq(lam0+1,lam):
            out=out*1j*K*get_a(lamt,mu,d)/jp[lamt]
        return out
    elif lam0==1 and lam==0 and mu==0:
        out=1j/(get_a(1,0,d)*K)
        return out
    elif lam0==2 and lam==0 and mu==0:
        out=W[lam0]*(-get_a(2,0,d)/get_a(1,0,d))
        return out
    elif lam0>2 and lam==0 and mu==0:
        out=W[lam0]*(-get_a(2,0,d)/get_a(1,0,d))
        for lamt in seq(lam+3,lam0):
            out=out*1j*K*get_a(lamt,mu,d)/jm[lamt-1]
        return out       
    elif lam0>lam and lam>0:
        out=W[lam0]
        for lamt in seq(lam+1,lam0):
            out=out*1j*K*get_a(lamt,mu,d)/jm[lamt-1]
        return out


# $$
# \lim_{p\to0}\partial_{p}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\begin{cases}
# \partial_{p}W_{\lambda_{0}}^{\mu} & \lambda_{0}=\lambda\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}\frac{\partial_{p}w_{\tilde{\lambda}}^{\mu\left(+\right)}}{w_{\tilde{\lambda}}^{\mu\left(+\right)}}\right)\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\right)} & \lambda>\lambda_{0}\\
# \frac{-i\left(2+\left(a_{\lambda_{0}+1}^{\mu}K\right)^{2}w_{\lambda_{0}+1}^{\mu\left(+\right)}\right)}{\left(a_{\lambda_{0}}^{\mu}K\right)^{3}} & \lambda_{0}=1,\ \lambda=\mu=0\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}-W_{\lambda_{0}}^{\mu}\frac{2}{\left(a_{1}^{0}K\right)^{2}}\right)\left(-\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\right) & \lambda_{0}=2,\ \lambda=\mu=0\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\left(\sum_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}-\frac{2}{\left(a_{1}^{0}K\right)^{2}}\right)\right)\left(-\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\right)\left(\prod_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}\right) & \lambda_{0}>2,\ \lambda=\mu=0\\
# \left(W_{\lambda_{0}}^{\mu}\partial_{p}w_{1}^{\mu\left(-\right)}\right)iKa_{2}^{\mu} & \lambda_{0}=2,\ \lambda=1,\ \mu=0\\
# \left(W_{\lambda_{0}}^{\mu}\partial_{p}w_{1}^{\mu\left(-\right)}\right)iKa_{2}^{\mu}\left(\prod_{\tilde{\lambda}=\lambda+2}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}\right) & \lambda_{0}>2,\ \lambda=1,\mu=0\\
# \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\right)}}\right)\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\right)} & \lambda_{0}>\lambda>1
# \end{cases}
# $$

# In[38]:

@jit(nopython=True)
def get_dG_zero(lam,lam0,mu,K,jm,jp,W,dW,dwm,dwp,d=3):
    if lam==lam0:
        return dW[lam0]
    elif lam0<lam:
        Sum=0.0
        Prod=1.0
        for lamt in seq(lam0+1,lam):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jp[lamt]
            Sum=Sum+dwp[lamt] * jp[lamt]
        return ((dW[lam0] + W[lam0]*Sum)*Prod)
    elif lam0==1 and lam==0 and mu==0:
        out=-1j*(2.0+(get_a(lam0+1,mu,d)*K)**2 *(1.0/jp[lam0+1]))/((get_a(lam0,mu,d)*K)**3)
        return out
    elif lam0==2 and lam==0 and mu==0:
        out=(dW[lam0]-2.0*W[lam0]/((get_a(1,0,d)*K)**2))*(-get_a(2,0,d)/get_a(1,0,d))
        return out
    elif lam0>2 and lam==0 and mu==0:
        Sum=-2.0/((get_a(1,0,d)*K)**2)
        Prod=1.0
        for lamt in seq(lam+3,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jm[lamt-1]
            Sum=Sum+dwm[lamt-1] * jm[lamt-1]
        return ((dW[lam0] + W[lam0]*Sum)*Prod)
    elif lam0==2 and lam==1 and mu==0:
        return W[lam0]*dwm[1]*1j*K*get_a(2,0,d)
    elif lam0>2 and lam==1 and mu==0:
        Prod=W[lam0]*dwm[1]*1j*K*get_a(2,0,d)
        for lamt in seq(lam+2,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jm[lamt-1]   
        return Prod
    elif lam0>lam and lam>1:
        Sum=0.0
        Prod=1.0
        for lamt in seq(lam+1,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jm[lamt-1]
            Sum=Sum+dwm[lamt-1] * jm[lamt-1]
        return ((dW[lam0] + W[lam0]*Sum)*Prod)  


# ## Residue Evaluation

# We would like to calculate first order residues defined as
# \[
# Res\left(\epsilon_{l}\right)=\lim_{p\to\epsilon_{l}}\left[\left(p-\epsilon_{l}\right)\mathcal{G}_{l}\right]
# \]
# by L'Hospital's Rule
# \[
# Res\left(\epsilon_{l}\right)=\frac{1}{\left.\frac{\partial}{\partial p}\left(\frac{1}{\mathcal{G}\left(p\right)}\right)\right|_{p=\epsilon_{l}}}
# \]

# In[134]:

# Calculate Residues using continued fraction
# Inputs:
#     K = wavelength (positive real number)
#     eig = eigenvalue at which to evaluate (complex scalar)
#     mu = z-component (intiger)
#     nlam = number of lambda values to return
#     lamMax = accurady, must be > nlam
# Outputs:
#      Res = Res_lam0,lam (nlam x nlam numpy array)
@jit(nopython=True)
def largeKResidues(K,eig,mu,nlam=10,lamMax=500,d=3):
    p=eig
    jp=get_jp(p,mu,K,lamMax)
    jm=get_jm(p,mu,K,lamMax)
    djp=get_djp(p,mu,K,lamMax,jp)
    djm=get_djm(p,mu,K,lamMax,jm)
    #W=get_W(p,mu,K,lamMax,jp,jm)
    Res=np.zeros((nlam,nlam),dtype=type(1+1j))*np.NaN
    for lam in range(mu,nlam):
        for lam0 in range(mu,nlam):
            Res[lam0,lam]=1.0/InvGprimeGeneral(p,lam,lam0,mu,K,jp,jm,djp,djm,d)
            #IG=InvG(p,lam,lam0,mu,K,jp,jm,d=3)
            #tol=0.0001
            #if (abs(IG)>tol):
            #    print('Error: Pole at p=',p,' 1/G=',IG,' > tol=',tol)
            #    print('lam',lam,'lam0',lam0,'mu',mu)
            #    return np.NaN
            #check=InvGprime(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3)
            #if (abs(check-Res[lam0,lam])>tol):
            #    print('Error: Pole at p=',p,' methode1',Res[lam0,lam],\
            #          ' methode2',check)
            #    return np.NaN                
    return Res


# \[
# C_{\lambda,l}^{\mu}=\begin{cases}
# 1 & \lambda=l\\
# \left(-iK\right)^{l-\lambda}\prod_{j=\lambda}^{l-1}\frac{a_{j+1}^{\mu}}{l\left(l+1\right)-j\left(j+1\right)} & l>\lambda\\
# \left(-iK\right)^{\lambda-l}\prod_{j=l+1}^{\lambda}\frac{a_{j}^{\mu}}{l\left(l+1\right)-j\left(j+1\right)} & l<\lambda
# \end{cases}
# \]

# In[135]:

@jit(nopython=True)
def get_C(K,l,lam,mu,d=3):
    prod=1
    Wi=l*(l+d-2)
    if l>lam:
        for j in seq(lam,l-1):  # lam to l-1
            Wj=j*(j+d-2)
            prod=prod*get_a(j+1,mu,d)/(Wi-Wj)
        prod=prod*(-1j*K)**(l-lam)
    elif l<lam:
        for j in seq(l+1,lam): #  l+1 to lam
            Wj=j*(j+d-2)
            prod=prod*get_a(j,mu,d)/(Wi-Wj)
        prod=prod*(-1j*K)**(lam-l)
    return prod


# \[
# \lim_{p\to\epsilon_{l}}\left[\left(p-\epsilon_{l}\right)\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(K,p\right)\right]\approx C_{\lambda_{0},l}^{\mu}C_{\lambda,l}^{\mu}
# \]

# In[136]:

# l = Eigenvalue number
# lam1, lam2: spherical harmonic l index
# mu: spherical harmonic mu index
@jit(nopython=True)
def SmallAysmpRes(K,l,lam1,lam2,mu,d=3):
    Res1=get_C(K,l,lam1,mu,d)
    Res2=get_C(K,l,lam2,mu,d)
    return Res1*Res2


# In[137]:

@jit(nopython=True)
def residues(K,eig,l,mu,nlam=10,d=3,lamMax=500,cutoff=10**-11):
    # use small K limit if K<10^-3 or the limit < cutoff
    if K<10**-3:
        res=np.zeros((nlam,nlam),dtype=type(1+1j))*np.NaN
        for lam in range(mu,nlam):
            for lam0 in range(mu,nlam):
                #print('lam0=',lam0,'lam',lam,'mu',mu)
                smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
                res[lam0,lam]=smallK
        return res   
    
    res=largeKResidues(K,eig,mu,nlam,lamMax=lamMax)
    for lam in range(mu,nlam):
        for lam0 in range(mu,nlam):
            #print('lam0=',lam0,'lam',lam,'mu',mu)
            smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
            if abs(smallK)<cutoff:
                res[lam0,lam]=smallK
    return res


# In[ ]:

@jit(nopython=True)
def bothResidues(K,eig,l,mu,nlam=10,d=3):
    cutoff=10**-11
    lamMax=500
    res=largeKResidues(K,eig,mu,nlam,lamMax=lamMax)
    resSmall=res*np.NaN
    for lam in range(mu,nlam):
        for lam0 in range(mu,nlam):
            #print('lam0=',lam0,'lam',lam,'mu',mu)
            resSmall[lam0,lam]=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
    return res, resSmall

