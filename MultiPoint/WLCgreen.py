

# ============================================================ #
#
# # Eigenvalues and Summand for WLC
#
# ============================================================ #


"""
This is a python implementation of eigenvalue and summand evaluation for the worm like chain propagator based on the 2008 PRE paper "End-to-End distribution for a wormlike chain in arbitrary dimensions" by Shafigh Mehraeen, Bariz Sudhanshu, Elena Koslover, and Andrew Spakowitz.
In also contains code for calculating higher derivatives of the poles.
This code was extended and translated into python from matlab by Quinn MacPherson in 2017.
"""



import numpy as np
from MultiPoint.utils import sphinx_compat_jit as jit
#from numba import jit




@jit(nopython=True)
def get_a(lam,mu,d=3):
    """
    coding: utf-8

    .. math::
        a_j^\mu=\sqrt{\\frac{(j-\mu)(j+\mu+D-3)}{(2j+D-2)(2j+D-4)}}

    """
    return np.sqrt((lam-mu)*(lam+mu+d-3)/ float((2*lam+d-2)*(2*lam+d-4)))

get_a(1,1,3)

# --------------------------------------------------- #
# ## Eigenvalues
# --------------------------------------------------- #

# ### Intermediate K
#_____________________________#




def IntermediateKEigenValues(K,ORDEig,mu,d=3):
    """
    Returns top ORDEig eigenvalus using Intermediate K methode.

    

    We would like to calculate the eigenvalues denoted :math:`\\epsilon_l` in the paper. These values depend on the z-component quantum number :math:`\mu` and the magnitude of the Fourier frequency :math:`K`.

    In the intermediate K regime we use the method in appendix A.2. We find the eigenvalues of the tri-diagonal matrix :math:`J` where the diagonal elements are 

    .. math::
        J_{j,j}=(\mu+j-1)(\mu+j+D-3)

    And off diagonal elements by

    .. math::
        J_{j+1,j}=J_{j,j+1}=-ia_{\mu+j}^{\mu}K

    only the eigenvalues ORDEig eigenvalues with the largest real part are returned.

    """
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
    out[abs(mu):]=evals[0:ORDEig-abs(mu)]
    return out
    



@jit(nopython=True)
def make_E(n_rows,mu,K,d=3):
    """
    E matrix for intermediate K. 

    """
    E=np.zeros((n_rows,n_rows),dtype=type(1+1j))
    left=(0+0j)*np.NaN
    for row in range(0,n_rows):
        j=row+1
        if K<=1:
            right=-1j*K*get_a(j+mu,mu,d) # for above diagonal
            if row>0:
                left=-1j*K*get_a(j+mu-1,mu,d) # for below diagonal
            diag=(mu+j-1)*(mu+j+d-3) # diagonal element
            if row==0:
                E[row,0:2]=[diag,right]
            elif row==n_rows-1:
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
            elif row==n_rows-1:
                E[row,row-1:row+1]=[left,diag]
            else:
                E[row,row-1:row+2]=[left,diag,right]
    return E






# ### Large K
#_____________________________#


def Epsilon(lam,d,alpha,mu):
    """
    implementation equivalent to table II for E from PRE 77, 061803 (2008)

    """
    beta=-np.sqrt(2)/4*(1+1j);
    m=mu+(d-3)/2.0; # Q.J.M. changes this line 8/1/15
    n=2*lam+m+1; # also know as s

    epsilon_0=(-1.0/2/beta)**(-1)*(n/2.0);
    epsilon_1=(-1.0/2/beta)**( 0)*(-1.0/8*(n**2+3-3*m**2)-m*(m+1));
    epsilon_2=(-1.0/2/beta)**( 1)*(-1.0/2**5*n*(n**2+3-9*m**2));
    epsilon_3=(-1.0/2/beta)**( 2)*(-1.0/2**8*(5*n**4+34*n**2+9)-(102*n**2+42)*m**2+33*m**4);
    epsilon_4=(-1.0/2/beta)**( 3)*(-1.0/2**11*n*(33*n**4+410*n**2+405)-(1230*n**2+1722)*m**2 +813*m**4);
    epsilon_5=(-1.0/2/beta)**( 4)*(-1.0/2**12*9*(7*n**6+140*n**4+327*n**2 +54-(420*n**4+1350*n**2+286)*m**2 +(495*n**2+314)*m**4-82*m**6));

    value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha**2+ epsilon_4*alpha**3+epsilon_5*alpha**4;
    return value




def LargeKEigenValues(k,ORDEig,mu,d=3):
    """
    Large K eigenvalues are implemented via

    .. math::
        \\epsilon_{2l}=iK-\mu(\mu+D-2)-\\frac{1}{\\alpha}\sum_{n=0}^{\infty}E_n^{(l)}\\alpha^n

    """
    Eig=np.zeros(ORDEig+1,dtype='complex')*np.NaN
    for l in range(0,ORDEig,2):
        alpha=1.0/np.sqrt(8*k)
        Eig[l]=1j*k-mu*(mu+d-2)-Epsilon(l/2.0,d,alpha,mu)
        Eig[l+1]=np.conj(Eig[l])
        
    out=np.zeros(ORDEig,'complex')*np.NaN
    out[mu:]=Eig[0:ORDEig-mu]
    return out
    



def wlc_eigvals(k,ORDEig,mu,d=3,cutoff=800):
    """
    Return eigenvalues calculated by Intermediate or Large K procedure

    Cutoff based on graph

    """
    if abs(k)<cutoff:
        Eig=IntermediateKEigenValues(k,ORDEig,mu,d)
    else:
        Eig=LargeKEigenValues(k,ORDEig,mu,d)
    return Eig


# --------------------------------------------------- #
# ## GLM
# --------------------------------------------------- #

# ### G
#_____________________________#


@jit(nopython=True)
def seq(first,last):
    """
    sequence of integers from first to last

    for example seq(2,5)=2,3,4,5

    """
    if first<=last:
        return np.arange(first,last+1)
    else:
        return np.arange(first,last-1,-1)




@jit(nopython=True)
def get_P(lam,p,d=3):
    """
    .. math::
        P_{\lambda}=p+\lambda\left(\lambda+d-2\\right)

    """
    return p+lam*(lam+d-2)


# #### j
#_____________________________#



@jit(nopython=True)
def get_jp(p,mu,K,lamMax,d=3):
    """
    Calculate j(+), output indexed by lambda size: lamMax+1

    .. math::
        j_{\lambda}^{\mu\left(+\\right)}=P_{\lambda}+\left(a_{\lambda+1}^{\mu}K\\right)^{2}\\frac{1}{j_{\lambda+1}^{\mu\left(+\\right)}}

    """
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




@jit(nopython=True)
def get_jm(p,mu,K,lamMax,d=3):
    """
    .. math::
        j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        P_{\lambda}+\left(a_{\lambda}^{\mu}K\\right)^{2}\\frac{1}{j_{\lambda-1}^{\mu\left(-\\right)}}
        &\ \ \ & \lambda\geq\mu+1 \\\\
        P_{\lambda} &\ \ \ & \lambda=\mu

        \\end{cases}

    """
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




@jit(nopython=True)
def get_djp(p,mu,K,lamMax,jp,d=3):
    """
    .. math::
        \\frac{\partial}{\partial p}j_{\lambda}^{\mu\left(+\\right)}=1-\left(a_{\lambda+1}^{\mu}K\\right)^{2}\\frac{1}{\left(j_{\lambda+1}^{\mu\left(+\\right)}\\right)^{2}}\\frac{\partial}{\partial p}j_{\lambda+1}^{\mu\left(+\\right)}

    """
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




@jit(nopython=True)
def get_djm(p,mu,K,lamMax,jm,d=3):
    """
    .. math::
        \\frac{\partial}{\partial p}j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        1-\left(a_{\lambda}^{\mu}K\\right)^{2}\\frac{1}{\left(j_{\lambda-1}^{\mu\left(-\\right)}\\right)^{2}}\\frac{\partial}{\partial
        p}j_{\lambda-1}^{\mu\left(-\\right)} &\ \ \ & \lambda\geq\mu+1 \\\\
        1 &\ \ \ & \lambda=\mu

        \\end{cases}

    """
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




@jit(nopython=True)
def get_ddjp(p,mu,K,lamMax,jp,djp,d=3):
    """
    .. math::
        d^{2}j_{\lambda}^{\mu\left(+\\right)}=\left(a_{\lambda+1}^{\mu}K\\right)^{2}\left(\\frac{2\left(dj_{\lambda+1}^{\mu\left(+\\right)}\\right)^{2}-\left(d^{2}j_{\lambda+1}^{\mu\left(+\\right)}\\right)j_{\lambda+1}^{\mu\left(+\\right)}}{\left(j_{\lambda+1}^{\mu\left(+\\right)}\\right)^{3}}\\right)

    """
    ddjp=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    ddjp[lamMax]= 1.0 # Just guessing
    for lam in seq(lamMax-1,mu):
        j=jp[lam+1]
        dj=djp[lam+1]
        ddj=ddjp[lam+1]
        ddjp[lam] = (get_a(lam+1,mu,d=d)*K)**2 * (2*(dj**2)-ddj*j)/(j**3)
    return ddjp 




@jit(nopython=True)
def get_ddjm(p,mu,K,lamMax,jm,djm,d=3):
    """
    .. math::
        d^{2}j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        \left(a_{\lambda}^{\mu}K\\right)^{2}\left(\\frac{2\left(dj_{\lambda-1}^{\mu\left(-\\right)}\\right)^{2}-\left(d^{2}j_{\lambda-1}^{\mu\left(-\\right)}\\right)j_{\lambda-1}^{\mu\left(-\\right)}}{\left(j_{\lambda-1}^{\mu\left(-\\right)}\\right)^{3}}\\right)
        &\ \ \ & \lambda\geq\mu+1 \\\\
        0 &\ \ \ & \lambda=\mu

        \\end{cases}

    """
    ddjm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    ddjm[mu]= 0.0
    for lam in seq(mu+1,lamMax):
        j=jm[lam-1]
        dj=djm[lam-1]
        ddj=ddjm[lam-1]
        ddjm[lam] =(get_a(lam,mu,d=d)*K)**2 * (2*(dj**2)-ddj*j)/(j**3)
    return ddjm




@jit(nopython=True)
def get_dddjp(p,mu,K,lamMax,jp,djp,ddjp,d=3):
    """
    .. math::
        d^{3}j_{\lambda}^{\mu\left(+\\right)}=\left(a_{\lambda+1}^{\mu}K\\right)^{2}\left(\\frac{6\left(d^{2}j_{\lambda+1}^{\mu\left(+\\right)}\\right)\left(dj_{\lambda+1}^{\mu\left(+\\right)}\\right)j_{\lambda+1}^{\mu\left(+\\right)}-6\left(dj_{\lambda+1}^{\mu\left(+\\right)}\\right)^{3}-\left(d^{3}j_{\lambda+1}^{\mu\left(+\\right)}\\right)\left(j_{\lambda+1}^{\mu\left(+\\right)}\\right)^{2}}{\left(j_{\lambda+1}^{\mu\left(+\\right)}\\right)^{4}}\\right)

    """
    dddjp=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    dddjp[lamMax]= 1.0 # Just guessing
    for lam in seq(lamMax-1,mu):
        j=jp[lam+1]
        dj=djp[lam+1]
        ddj=ddjp[lam+1]
        dddj=dddjp[lam+1]
        dddjp[lam] = (get_a(lam+1,mu,d=d)*K)**2 * (6.0*ddj*dj*j-6.0*(dj**3)-dddj*(j**2))/(j**4)
    return dddjp 




@jit(nopython=True)
def get_dddjm(p,mu,K,lamMax,jm,djm,ddjm,d=3):
    """
    .. math::
        d^{3}j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        \left(a_{\lambda}^{\mu}K\\right)^{2}\left(\\frac{6\left(d^{2}j_{\lambda-1}^{\mu\left(-\\right)}\\right)\left(dj_{\lambda-1}^{\mu\left(-\\right)}\\right)j_{\lambda-1}^{\mu\left(-\\right)}-6\left(dj_{\lambda-1}^{\mu\left(-\\right)}\\right)^{3}-\left(d^{3}j_{\lambda-1}^{\mu\left(-\\right)}\\right)\left(j_{\lambda-1}^{\mu\left(-\\right)}\\right)^{2}}{\left(j_{\lambda-1}^{\mu\left(-\\right)}\\right)^{4}}\\right)
        &\ \ \ & \lambda\geq\mu+1 \\\\
        0 &\ \ \ & \lambda=\mu

        \\end{cases}

    """
    dddjm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    dddjm[mu]= 0.0
    for lam in seq(mu+1,lamMax):
        j=jm[lam-1]
        dj=djm[lam-1]
        ddj=ddjm[lam-1]
        dddj=ddjm[lam-1]
        dddjm[lam] =(get_a(lam,mu,d=d)*K)**2 * (6.0*ddj*dj*j-6.0*(dj**3)-dddj*(j**2))/(j**4)
    return dddjm


# #### w
#_____________________________#


@jit(nopython=True)
def get_w(j,mu):
    """
    Calculate :math:`w=\\frac{1}{j}`

    """
    j[j==0.0+0.0j]=np.NaN
    w=1.0/j
    return w




@jit(nopython=True)
def get_dw(j,dj,mu):
    """
    .. math::
        \partial_{p}w_{\lambda}^{\mu\left(\pm\\right)}=\\frac{-1}{\left(j_{\lambda}^{\mu\left(\pm\\right)}\\right)^{2}}\partial_{p}j_{\lambda}^{\mu\left(\pm\\right)}

    """
    j[j==0.0+0.0j]=np.NaN
    dw=-dj/(j**2)
    #j[j==0.0+0.0j]=np.NaN
    #dw=j*np.NaN
    #dw[mu:]=-dj[mu:]/(j[mu:]**2)
    return dw




@jit(nopython=True)
def get_ddw(j,dj,ddj,mu):
    """
    .. math::
        \\frac{\partial^{2}}{\partial p^{2}}w_{\lambda}^{\mu\left(\pm\\right)}=\left(\\frac{2\left(dj_{\lambda}^{\mu\left(\pm\\right)}\\right)^{2}-\left(d^{2}j_{\lambda}^{\mu\left(\pm\\right)}\\right)j_{\lambda}^{\mu\left(\pm\\right)}}{\left(j_{\lambda}^{\mu\left(\pm\\right)}\\right)^{3}}\\right)

    """
    #dw=-dj/(j**2)
    j[j==0.0+0.0j]=np.NaN
    ddw=j*np.NaN
    ddw[mu:]=(2*(dj[mu:]**2)-ddj[mu:]*j[mu:])/(j[mu:]**3)
    return ddw




@jit(nopython=True)
def get_dddw(j,dj,ddj,dddj,mu):
    """
    .. math::
        \\frac{\partial^{3}}{\partial p^{3}}w_{\lambda}^{\mu\left(\pm\\right)}=\left(\\frac{6\left(d^{2}j_{\lambda}^{\mu\left(\pm\\right)}\\right)\left(dj_{\lambda}^{\mu\left(\pm\\right)}\\right)j_{\lambda}^{\mu\left(\pm\\right)}-6\left(dj_{\lambda}^{\mu\left(\pm\\right)}\\right)^{3}-d^{3}j_{\lambda}^{\mu\left(\pm\\right)}\left(j_{\lambda}^{\mu\left(\pm\\right)}\\right)^{2}}{\left(j_{\lambda}^{\mu\left(\pm\\right)}\\right)^{4}}\\right)

    """
    #dw=-dj/(j**2)
    j[j==0.0+0.0j]=np.NaN
    dddw=j*np.NaN
    dddw[mu:]= (6.0*ddj[mu:]*dj[mu:]*j[mu:]-6.0*(dj[mu:]**3)-dddj[mu:]*(j[mu:]**2))/(j[mu:]**4)
    return dddw




@jit(nopython=True)
def get_W(p,mu,K,lamMax,jp,jm,d=3):
    """
    .. math::
        W_{\lambda}^{\mu}=\\frac{1}{\left(a_{\lambda}^{\mu}K\\right)^{2}w_{\lambda-1}^{\mu\left(-\\right)}+P+\left(a_{\lambda+1}^{\mu}K\\right)^{2}w_{\lambda+1}^{\mu\left(+\\right)}}

    """
    W=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    lam=mu
    W[lam]=get_P(lam,p,d=d)+(get_a(lam+1,mu,d=d)*K)**2 / jp[lam+1]
    for lam in seq(mu+1,lamMax-1):
        W[lam]=(get_a(lam,mu,d=d)*K)**2 /jm[lam-1] + get_P(lam,p,d=d)+ (get_a(lam+1,mu,d=d)*K)**2 / jp[lam+1]
    W[W==0.0+0.0j]=np.NaN
    return (1.0/W)




@jit(nopython=True)
def get_dW(p,mu,K,lamMax,dwp,dwm,W,d=3):
    """
    .. math::
        \partial_{p}W_{\lambda}^{\mu}=-\left(W_{\lambda}^{\mu}\\right)^{2}\left(\left(a_{\lambda}^{\mu}K\\right)^{2}\partial_{p}w_{\lambda-1}^{\mu\left(-\\right)}+1+\left(a_{\lambda+1}^{\mu}K\\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\\right)}\\right)

    """
    dW=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    lam=mu
    dW[lam]=1.0 + (get_a(lam+1,mu,d=d)*K)**2 * dwp[lam+1]
    for lam in seq(mu+1,lamMax-1):
        dW[lam]=(get_a(lam,mu,d=d)*K)**2 *dwm[lam-1] + 1.0 + (get_a(lam+1,mu,d=d)*K)**2 * dwp[lam+1]
    return (-dW*(W**2))


# #### G
#_____________________________#



@jit(nopython=True)
def get_G(p,lam,lam0,mu,K,jm,jp,W,d=3):
    """
    .. math::
        \mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\\begin{cases}

        W_{\lambda_{0}}^{\mu} &\ \ \ & \lambda_{0}=\lambda \\\\
        W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\\right)}
        && \lambda_{0}<\lambda \\\\
        W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}
        && \lambda_{0}>\lambda

        \\end{cases}

    """
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




@jit(nopython=True)
def get_dG(p,lam,lam0,mu,K,jm,jp,W,dW,dwm,dwp,d=3):
    """
    .. math::
        \partial_{p}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\\begin{cases}

        \partial_{p}W_{\lambda_{0}}^{\mu} &\ \ \ & \lambda_{0}=\lambda \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}\\frac{\partial_{p}w_{\tilde{\lambda}}^{\mu\left(+\\right)}}{w_{\tilde{\lambda}}^{\mu\left(+\\right)}}\\right)\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\\right)}
        && \lambda_{0}<\lambda \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}\\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}\\right)\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}
        && \lambda_{0}>\lambda

        \\end{cases}

    """
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
#_____________________________#

# #### X1
#_____________________________#



@jit(nopython=True)
def get_X1(p,lam0,mu,K,jp,jm,d=3):
    """
    .. math::
        X_{1}=\left(a_{\lambda_{0}}^{\mu}K\\right)^{2}w_{\lambda_{0}-1}^{\mu\left(-\\right)}+P_{\lambda_{0}}+\left(a_{\lambda_{0}+1}^{\mu}K\\right)^{2}w_{\lambda_{0}+1}^{\mu\left(+\\right)}

    """
    if lam0 == mu:
        out= 0 + get_P(lam0,p,d=d) +(get_a(lam0+1,mu,d=d)*K)**2 * 1/jp[lam0+1]
    else:
        out= (get_a(lam0,mu,d=d)*K)**2 * 1/jm[lam0-1] + get_P(lam0,p,d=d) +(get_a(lam0+1,mu,d=d)*K)**2 * 1/jp[lam0+1]
    return out




@jit(nopython=True)
def get_dX1(lam0,mu,K,dwp,dwm,d=3):
    """
    .. math::
        \\frac{\partial}{\partial p}X_{1}=\left(a_{\lambda_{0}}^{\mu}K\\right)^{2}\\frac{\partial}{\partial p}w_{\lambda_{0}-1}^{\mu\left(-\\right)}+1+\left(a_{\lambda_{0}+1}^{\mu}K\\right)^{2}\\frac{\partial}{\partial p}w_{\lambda_{0}+1}^{\mu\left(+\\right)}

    """
    if lam0 == mu:
        dX1= 0 + 1 +(get_a(lam0+1,mu,d=d)*K)**2 * dwp[lam0+1] 
    else:
        dX1= (get_a(lam0,mu,d=d)*K)**2 * dwm[lam0-1] + 1 +(get_a(lam0+1,mu,d=d)*K)**2 * dwp[lam0+1] 
    return dX1




@jit(nopython=True)
def get_ddX1(lam0,mu,K,ddwp,ddwm,d=3):
    """
    .. math::
        \\frac{\partial^{2}}{\partial p^{2}}X_{1}=\left(a_{\lambda_{0}}^{\mu}K\\right)^{2}\\frac{\partial^{2}}{\partial p^{2}}w_{\lambda_{0}}^{\mu\left(-\\right)}+\left(a_{\lambda_{0}+1}^{\mu}K\\right)^{2}\\frac{\partial^{2}}{\partial p^{2}}w_{\lambda_{0}+1}^{\mu\left(+\\right)}

    """
    if lam0 == mu:
        ddX1= 0 + 0 +(get_a(lam0+1,mu,d=d)*K)**2 * ddwp[lam0+1] 
    else:
        ddX1= (get_a(lam0,mu,d=d)*K)**2 * ddwm[lam0-1] + 0 +(get_a(lam0+1,mu,d=d)*K)**2 * ddwp[lam0+1] 
    return ddX1




@jit(nopython=True)
def get_dddX1(lam0,mu,K,dddwp,dddwm,d=3):
    """
    .. math::
        \\frac{\partial^{3}}{\partial p^{3}}X_{1}=\left(a_{\lambda_{0}}^{\mu}K\\right)^{2}\\frac{\partial^{3}}{\partial p^{3}}w_{\lambda_{0}}^{\mu\left(-\\right)}+\left(a_{\lambda_{0}+1}^{\mu}K\\right)^{2}\\frac{\partial^{3}}{\partial p^{3}}w_{\lambda_{0}+1}^{\mu\left(+\\right)}

    """
    if lam0 == mu:
        dddX1= 0 + 0 +(get_a(lam0+1,mu,d=d)*K)**2 * dddwp[lam0+1] 
    else:
        dddX1= (get_a(lam0,mu,d=d)*K)**2 * dddwm[lam0-1] + 0 +(get_a(lam0+1,mu,d=d)*K)**2 * dddwp[lam0+1] 
    return dddX1


# #### Xpm
#_____________________________#



@jit(nopython=True)
def get_Sp(lam,lam0,jp,djp,pm):
    """
    .. math::
        S^{\left(\pm\\right)}=\sum_{n=1}^{\left|\lambda-\lambda_{0}\\right|}\\frac{\partial_{p}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}{j_{\lambda_{0}\pm n}^{\mu\left(+\\right)}}

    """
    S=0.0+0.0j
    for n in seq(1,abs(lam-lam0)):
        j=jp[lam0+pm*n]
        dj=djp[lam0+pm*n]
        S=S+(dj/j)
    return S




@jit(nopython=True)
def get_dSp(lam,lam0,jp,djp,ddjp,pm):
    """
    .. math::
        \\frac{\partial}{\partial p}S^{\left(\pm\\right)}=\sum_{n=1}^{\left|\lambda-\lambda_{0}\\right|}\left(\\frac{\partial_{p}^{2}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}{j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}-\\frac{\left(\partial_{p}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\\right)^{2}}{\left(j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\\right)^{2}}\\right)

    """
    S=0.0+0.0j
    for n in seq(1,abs(lam-lam0)):
        j=jp[lam0+pm*n]
        dj=djp[lam0+pm*n]
        ddj=ddjp[lam0+pm*n]
        S=S+((ddj/j)-(dj**2)/(j**2))
    return S




@jit(nopython=True)
def get_ddSp(lam,lam0,jp,djp,ddjp,dddjp,pm):
    """
    .. math::
        \\frac{\partial^{2}}{\partial p^{2}}S^{\left(+\\right)}=\sum_{n=1}^{\left|\lambda-\lambda_{0}\\right|}\left(\\frac{\partial_{p}^{3}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}{j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}+\\frac{2\left(\partial_{p}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\\right)^{3}}{\left(j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\\right)^{3}}-\\frac{3\partial_{p}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\partial_{p}^{2}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}}{\left(j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}\\right)^{2}}\\right)

    """
    S=0.0+0.0j
    for n in seq(1,abs(lam-lam0)):
        j=jp[lam0+pm*n]
        dj=djp[lam0+pm*n]
        ddj=ddjp[lam0+pm*n]
        dddj=dddjp[lam0+pm*n]
        S=S+((dddj/j)+2*(dj**3)/(j**3)-3*dj*ddj/(j**2))
    return S




@jit(nopython=True)
def get_Xp(lam,lam0,jp,pm):
    """
    .. math::
        X^{\left(\pm\\right)}=\prod_{n=1}^{\left|\lambda-\lambda_{0}\\right|}j_{\lambda_{0}\pm n}^{\mu\left(\pm\\right)}

    """
    prod=1.0+0.0j
    for n in seq(1,abs(lam-lam0)):
        prod=prod*jp[lam0+pm*n]
    return prod




@jit(nopython=True)
def get_dXp(S,X):
    """
    .. math::
        \\frac{\partial}{\partial p}X^{\left(\pm\\right)}=S^{\left(\pm\\right)}X^{\left(\pm\\right)}

    """
    return X*S




@jit(nopython=True)
def get_ddXp(S,dS,X):
    """
    .. math::
        \\frac{\partial^{2}}{\partial p^{2}}X^{\left(\pm\\right)}=\left(\partial_{p}S^{\left(\pm\\right)}\\right)X^{\left(\pm\\right)}+\left(S^{\left(\pm\\right)}\\right)^{2}X^{\left(\pm\\right)}

    """
    return (dS+(S**2))*X




@jit(nopython=True)
def get_dddXp(S,dS,ddS,X):
    """
    .. math::
        \\frac{\partial^{3}}{\partial p^{3}}X^{\left(\pm\\right)}=\left(\partial_{p}^{2}S^{\left(\pm\\right)}\\right)X^{\left(\pm\\right)}+3\left(\partial_{p}S^{\left(\pm\\right)}\\right)S^{\left(\pm\\right)}X^{\left(+\\right)}+\left(S^{\left(\pm\\right)}\\right)^{3}X^{\left(\pm\\right)}

    """
    return (ddS + 3*dS*S + S**3)*X




@jit(nopython=True)
def get_cp(lam,lam0,mu,K,d=3):
    """
    .. math::
        c^{\left(+\\right)}=\prod_{n=1}^{\lambda-\lambda_{0}}iKa_{\lambda_{0}+n}^{\mu}

    """
    prod=1
    for n in seq(1,lam-lam0):
        prod=prod*1j*K*get_a(lam0+n,mu,d)
    return prod




@jit(nopython=True)
def get_cm(lam,lam0,mu,K,d=3):
    """
    .. math::
        c^{\left(-\\right)}=\prod_{n=1}^{\lambda_{0}-\lambda}iKa_{\lambda_{0}+1-n}^{\mu}

    """
    prod=1
    for n in seq(1,lam0-lam):
        prod=prod*1j*K*get_a(lam0+1-n,mu,d)
    return prod


# #### 1/G
#_____________________________#



@jit(nopython=True)
def get_invG(cp,X1,Xp):
    """
    .. math::
        \\frac{1}{\mathcal{G}_{\lambda_{0}\lambda}^{\mu}\left(p\\right)}=\\begin{cases}

        X_{1} &\ \ \ & \lambda_{0}=\lambda \\\\
        \\frac{1}{c^{\left(+\\right)}}X_{1}X^{\left(+\\right)} && \lambda_{0}<\lambda \\\\
        \\frac{1}{c^{\left(-\\right)}}X_{1}X^{\left(-\\right)} && \lambda_{0}>\lambda

        \\end{cases}

    """
    return X1*Xp/cp




@jit(nopython=True)
def get_dinvG(cp,X1,dX1,Xp,dXp):
    """
    .. math::
        \\frac{\partial}{\partial p}\\frac{1}{\mathcal{G}_{\lambda_{0}\lambda}^{\mu}\left(p\\right)}=\\begin{cases}

        \partial_{p}X_{1} &\ \ \ & \lambda_{0}=\lambda \\\\
        \\frac{1}{c^{\left(+\\right)}}\left(X_{1}\partial_{p}X^{\left(+\\right)}+X^{\left(+\\right)}\partial_{p}X_{1}\\right)
        && \lambda_{0}<\lambda \\\\
        \\frac{1}{c^{\left(-\\right)}}\left(X_{1}\partial_{p}X^{\left(-\\right)}+X^{\left(-\\right)}\partial_{p}X_{1}\\right)
        && \lambda_{0}>\lambda

        \\end{cases}

    

    """
    return (X1*dXp+Xp*dX1)/cp




@jit(nopython=True)
def get_ddinvG(cp,X1,dX1,ddX1,Xp,dXp,ddXp):
    """
    .. math::
        \\frac{\partial^{2}}{\partial p^{2}}\\frac{1}{\mathcal{G}_{\lambda_{0}\lambda}^{\mu}\left(p\\right)}=\\begin{cases}

        \partial_{p}^{2}X_{1} &\ \ \ &  \lambda_{0}=\lambda \\\\
        \\frac{1}{c^{\left(+\\right)}}\left(X_{1}\partial_{p}^{2}X^{\left(+\\right)}+2\partial_{p}X_{1}\partial_{p}X^{\left(+\\right)}+X^{\left(+\\right)}\partial_{p}^{2}X_{1}\\right)
        && \lambda_{0}<\lambda \\\\
        \\frac{1}{c^{\left(-\\right)}}\left(X_{1}\partial_{p}^{2}X^{\left(-\\right)}+2\partial_{p}X_{1}\partial_{p}X^{\left(+\\right)}+X^{\left(-\\right)}\partial_{p}^{2}X_{1}\\right)
        && \lambda_{0}>\lambda

        \\end{cases}

    """
    return (X1*ddXp+2*dX1*dXp+Xp*ddX1)/cp




@jit(nopython=True)
def get_dddinvG(cp,X1,dX1,ddX1,dddX1, Xp,dXp,ddXp,dddXp):
    """
    .. math::
        \\frac{\partial^{3}}{\partial p^{3}}\\frac{1}{\mathcal{G}_{\lambda_{0}\lambda}^{\mu}\left(p\\right)}=\\begin{cases}

        \partial_{p}^{3}X_{1} &\ \ \ & \lambda_{0}=\lambda \\\\
        \\frac{1}{c^{\left(+\\right)}}\left(X_{1}\partial_{p}^{3}X^{\left(+\\right)}+3\partial_{p}X_{1}\partial_{p}^{2}X^{\left(+\\right)}+3\partial_{p}^{2}X_{1}\partial_{p}X^{\left(+\\right)}+X^{\left(+\\right)}\partial_{p}^{3}X_{1}\\right)
        && \lambda_{0}<\lambda \\\\
        \\frac{1}{c^{\left(-\\right)}}\left(X_{1}\partial_{p}^{3}X^{\left(-\\right)}+3\partial_{p}X_{1}\partial_{p}^{2}X^{\left(-\\right)}+3\partial_{p}^{2}X_{1}\partial_{p}X^{\left(-\\right)}+X^{\left(-\\right)}\partial_{p}^{3}X_{1}\\right)
        && \lambda_{0}>\lambda

        \\end{cases}

    """
    return (X1*dddXp + 3*dX1*ddXp + 3*ddX1*dXp + Xp*dddX1)/cp


# #### old
#_____________________________#


@jit(nopython=True)
def get_dw_old(j,dj,mu):
    """
    No longer in use.

    """
    dw=-dj/(j**2)
    #dw=j*np.NaN
    #dw[mu:]=-dj[mu:]/(j[mu:]**2)
    return dw



@jit(nopython=True)
def X1prime(p,lam0,mu,K,jp,jm,djp,djm,d=3):
    """
    Calculate :math:`X'_1`

    """
    if lam0 == mu:
        X1= 0 + 1 +(get_a(lam0+1,mu,d=d)*K)**2 * get_dw_old(jp[lam0+1],djp[lam0+1],mu) 
    else:
        X1= (get_a(lam0,mu,d=d)*K)**2 * get_dw_old(jm[lam0-1],djm[lam0-1],mu) + 1 +(get_a(lam0+1,mu,d=d)*K)**2 * get_dw_old(jp[lam0+1],djp[lam0+1],mu) 
    return X1



@jit(nopython=True)
def Xmprime(lam,lam0,jm,djm,Xm):
    """
    Calculate :math:`X'^{-}`
    """
    out=0.0
    for n in seq(1,lam0-lam):
        out=out+(djm[lam0-n]/jm[lam0-n])
    out=out*Xm
    return out



@jit(nopython=True)
def Xp_old(lam,lam0,jp):
    """
    No longer in use.

    """
    prod=1
    for n in seq(1,lam-lam0):
        prod=prod*jp[lam0+n]
    return prod



@jit(nopython=True)
def Xm_old(lam,lam0,jm):
    """
    No longer in use.

    """
    prod=1
    for n in seq(1,lam0-lam):
        prod=prod*jm[lam0-n]
    return prod



@jit(nopython=True)
def InvGprime(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3):
    """
    Calculate :math:`\\frac{1}{G'}`

    """
    if min(lam,lam0)<mu:
        return np.NaN
    tol=0.00001
    if abs(X1(p,lam0,mu,K,jp,jm)) > tol:
        print("Error in InvGpriem. X1 should be zero at pole")
        return np.NaN
    
    if lam0==lam:
        return X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)
    if lam0<lam:
        return (X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*Xp_old(lam,lam0,jp) /get_cp(lam,lam0,mu,K,d))
    if lam0>lam:
        return (X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*Xm_old(lam,lam0,jm) /get_cm(lam,lam0,mu,K,d))




@jit(nopython=True)
def InvGprimeGeneral(p,lam,lam0,mu,K,jp,jm,djp,djm,d=3):
    """
    .. math::
        \\frac{\partial}{\partial p}\left(\\frac{1}{\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(p\\right)}\\right)=\\begin{cases}

        \partial_{p}X_{1} &\ \ \ & \lambda_{0}=\lambda \\\\
        \\frac{X^{\left(+\\right)}\partial_{p}X_{1}}{c^{\left(+\\right)}} && \lambda_{0}<\lambda \\\\
        \\frac{X^{\left(-\\right)}\partial_{p}X_{1}}{c^{\left(-\\right)}} && \lambda_{0}>\lambda

        \\end{cases}

    because :math:`X_1=0` at the poles, :math:`\partial_p X^{ \\{\pm\\} }` does not generally diverge at the poles.

    """
    if min(lam,lam0)<mu:
        return np.NaN
    if lam0==lam:
        return X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)
    if lam0<lam:
        xp=Xp(lam,lam0,jp)
        return ((X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*xp +X1(p,lam0,mu,K,jp,jm)*Xpprime(lam,lam0,jp,djp,xp)) /get_cp(lam,lam0,mu,K,d))
    if lam0>lam:
        xm=Xm(lam,lam0,jm)
        return ((X1prime(p,lam,lam0,mu,K,jp,jm,djp,djm)*xm +X1(p,lam0,mu,K,jp,jm)*Xmprime(lam,lam0,jm,djm,xm)) /get_cm(lam,lam0,mu,K,d))


# ### Limit p -> 0
#_____________________________#



@jit(nopython=True)
def get_jm_zero(mu,K,lamMax,d=3):
    """
    .. math::
        \lim_{p\to0}j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        p &\ \ \ & \lambda=\mu=0 \\\\
        \left(a_{\lambda}^{\mu}K\\right)^{2}/p && \lambda=1,\ \mu=0 \\\\
        \lambda\left(\lambda+1\\right) && \lambda=2,\ \mu=0 \\\\
        \lambda\left(\lambda+1\\right) && \lambda=\mu>0 \\\\
        \lambda\left(\lambda+1\\right)+\left(a_{\lambda}^{\mu}K\\right)^{2}\\frac{1}{\lim_{p\to0}j_{\lambda-1}^{\mu\left(-\\right)}}
        && \mathrm{otherwise}

        \\end{cases}

    """
    jm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu!=0: 
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




@jit(nopython=True)
def get_wm_zero(jm,lamMax):
    """
    .. math::
        \lim_{p\to0}w_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        1/p &\ \ \  & \lambda=\mu=0 \\\\
        p/\left(a_{\lambda}^{\mu}K\\right)^{2} && \lambda=1,\ \mu=0 \\\\
        1/\lim_{p\to0}j_{\lambda}^{\mu\left(-\\right)} && \mathrm{\mathrm{otherwise}}

        \\end{cases}

    """
    wm=np.zeros(lamMax+1,dtype=type(1+1j))
    wm[0]=(1+0j)*np.Inf
    wm[2:lamMax+1]=1.0/jm[2:lamMax+1]




@jit(nopython=True)
def get_djm_zero(mu,K,lamMax,jm,d=3):
    """
    .. math::
        \lim_{p\to0}\\frac{\partial}{\partial p}j_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        1 &\ \ \ & \lambda=\mu=0 \\\\
        -\left(a_{\lambda}^{\mu}K\\right)^{2}/p^{2} && \lambda=1,\ \mu=0 \\\\
        1+\\frac{\left(a_{\lambda}^{\mu}\\right)^{2}}{\left(a_{\lambda-1}^{\mu}\\right)^{2}}
        && \lambda=2,\ \mu=0 \\\\
        1 && \lambda=\mu>0 \\\\
        1-\left(a_{\lambda}^{\mu}K\\right)^{2}\\frac{1}{j_{\lambda-1}^{\mu\left(-\\right)}{}^{2}}\\frac{\partial}{\partial
        p}j_{\lambda-1}^{\mu\left(-\\right)} && \mathrm{otherwise}

        \\end{cases}

    """
    djm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu!=0:
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




@jit(nopython=True)
def get_dwm_zero(jm,djm,K,lamMax,mu,d=3):
    """
    .. math::
        \lim_{p\to0}\partial_{p}w_{\lambda}^{\mu\left(-\\right)}=\\begin{cases}

        -\\frac{1}{p^{2}} &\ \ \ & \lambda=\mu=0 \\\\
        \\frac{1}{\left(a_{\lambda}^{\mu}K\\right)^{2}} && \lambda=1,\ \mu=0 \\\\
        \\frac{-1}{\left(j_{\lambda}^{\mu\left(-\\right)}\\right)^{2}}\partial_{p}j_{\lambda}^{\mu\left(-\\right)}
        && \mathrm{otherwise}

        \\end{cases}

    """
    dwm=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu!=0:
        dwm[mu:]=-1.0*djm[mu:]/(jm[mu:]**2)
    elif mu==0:
        dwm[0]=(-1.0+0j)*np.Inf
        dwm[1]=1/((get_a(1,0,d)*K)**2)
        dwm[2:lamMax+1]=-1.0*djm[2:lamMax+1]/(jm[2:lamMax+1]**2)
    return dwm




@jit(nopython=True)
def get_W_zero(mu,K,lamMax,jp,jm,d=3):
    """
    .. math::
        \lim_{p\to0}W_{\lambda}^{\mu}=\\begin{cases}

        \\frac{1}{\lambda\left(\lambda+1\\right)+\left(a_{\lambda+1}^{\mu}K\\right)^{2}w_{\lambda+1}^{\mu\left(+\\right)}}
        &\ \ \ & \lambda=0,2,\ \mu=0 \\\\
        \\frac{p}{\left(a_{\lambda}^{\mu}K\\right)^{2}} && \lambda=1,\ \mu=0 \\\\
        \\frac{1}{\left(a_{\lambda}^{\mu}K\\right)^{2}w_{\lambda-1}^{\mu\left(-\\right)}+\lambda\left(\lambda+1\\right)+\left(a_{\lambda+1}^{\mu}K\\right)^{2}w_{\lambda+1}^{\mu\left(+\\right)}}
        && \mathrm{otherwise}

        \\end{cases}

    """
    W=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu!=0:
        lam=mu
        W[lam]=get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1]
        for lam in seq(mu+1,lamMax-1):
            W[lam]=(get_a(lam,mu,d)*K)**2 /jm[lam-1] + get_P(lam,0,d)+ (get_a(lam+1,mu,d)*K)**2 / jp[lam+1]
        return (1/W)
    elif mu==0:
        lam=0
        W[0]=1.0/(get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        W[1]=0.0
        lam=2
        W[2]=1.0/(get_P(lam,0,d)+(get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        for lam in seq(3,lamMax-1):
            W[lam]=1.0/((get_a(lam,mu,d)*K)**2 /jm[lam-1] + get_P(lam,0,d)+ (get_a(lam+1,mu,d)*K)**2 / jp[lam+1])
        return W 




@jit(nopython=True)
def get_dW_zero(mu,K,lamMax,dwp,dwm,W,d=3):
    """
    .. math::
        \lim_{p\to0}\partial_{p}W_{\lambda}^{\mu}=\\begin{cases}

        -\left(W_{\lambda}^{\mu}\\right)^{2}\left(1+\left(a_{\lambda+1}^{\mu}K\\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\\right)}\\right)
        &\ \ \ & \lambda=0,\ \mu=0 \\\\
        \\frac{1}{\left(a_{\lambda}^{\mu}K\\right)^{2}} && \lambda=1,\ \mu=0 \\\\
        -\left(W_{\lambda}^{\mu}\\right)^{2}\left(\left(a_{\lambda}^{\mu}K\\right)^{2}\partial_{p}w_{\lambda-1}^{\mu\left(-\\right)}+1+\left(a_{\lambda+1}^{\mu}K\\right)^{2}\partial_{p}w_{\lambda+1}^{\mu\left(+\\right)}\\right)
        && \mathrm{otherwise}

        \\end{cases}

    """
    dW=np.zeros(lamMax+1,dtype=type(1+1j))*np.NaN
    if mu!=0:
        lam=mu
        dW[lam]=1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1]
        for lam in seq(mu+1,lamMax-1):
            dW[lam]=(get_a(lam,mu,d)*K)**2 *dwm[lam-1] + 1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1]
        return (-dW*(W**2))
    elif mu==0:
        lam=0
        dW[0]=-(W[lam]**2)*(1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1])
        dW[1]=(get_a(1,0,d)*K)**(-2)
        for lam in seq(2,lamMax-1):
            dW[lam]=-(W[lam]**2)*((get_a(lam,mu,d)*K)**2 *dwm[lam-1] + 1.0 + (get_a(lam+1,mu,d)*K)**2 * dwp[lam+1])
        return dW




@jit(nopython=True)
def get_G_zero(lam,lam0,mu,K,jm,jp,W,d=3):
    """
    .. math::
        \lim_{p\to0}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\\begin{cases}

        W_{\lambda_{0}}^{\mu} &\ \ \ & \lambda_{0}=\lambda \\\\
        W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\\right)}
        && \lambda_{0}<\lambda \\\\
        \\frac{i}{a_{1}^{0}K} && \lambda_{0}=1,\ \lambda=\mu=0 \\\\
        W_{\lambda_{0}}^{\mu}\left(-\\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\\right) && \lambda_{0}=2,\ \lambda=\mu=0 \\\\
        W_{\lambda_{0}}^{\mu}\left(-\\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\\right)\left(\prod_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}\\right)
        && \lambda_{0}>2,\ \lambda=\mu=0 \\\\
        W_{\lambda_{0}}^{\mu}\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}
        && \lambda_{0}>\lambda>0

        \\end{cases}

    """
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




@jit(nopython=True)
def get_dG_zero(lam,lam0,mu,K,jm,jp,W,dW,dwm,dwp,d=3):
    """
    .. math::
        \lim_{p\to0}\partial_{p}\mathcal{G}_{\lambda_{0}\lambda}^{\mu}=\\begin{cases}

        \partial_{p}W_{\lambda_{0}}^{\mu} &\ \ \ & \lambda_{0}=\lambda \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}\\frac{\partial_{p}w_{\tilde{\lambda}}^{\mu\left(+\\right)}}{w_{\tilde{\lambda}}^{\mu\left(+\\right)}}\\right)\prod_{\tilde{\lambda}=\lambda_{0}+1}^{\lambda}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}}^{\mu\left(+\\right)}
        && \lambda>\lambda_{0} \\\\
        \\frac{-i\left(2+\left(a_{\lambda_{0}+1}^{\mu}K\\right)^{2}w_{\lambda_{0}+1}^{\mu\left(+\\right)}\\right)}{\left(a_{\lambda_{0}}^{\mu}K\\right)^{3}}
        && \lambda_{0}=1,\ \lambda=\mu=0 \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}-W_{\lambda_{0}}^{\mu}\\frac{2}{\left(a_{1}^{0}K\\right)^{2}}\\right)\left(-\\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\\right)
        && \lambda_{0}=2,\ \lambda=\mu=0 \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\left(\sum_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}\\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}-\\frac{2}{\left(a_{1}^{0}K\\right)^{2}}\\right)\\right)\left(-\\frac{a_{2}^{\mu}}{a_{1}^{\mu}}\\right)\left(\prod_{\tilde{\lambda}=\lambda+3}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}\\right)
        && \lambda_{0}>2,\ \lambda=\mu=0 \\\\
        \left(W_{\lambda_{0}}^{\mu}\partial_{p}w_{1}^{\mu\left(-\\right)}\\right)iKa_{2}^{\mu}
        && \lambda_{0}=2,\ \lambda=1,\ \mu=0 \\\\
        \left(W_{\lambda_{0}}^{\mu}\partial_{p}w_{1}^{\mu\left(-\\right)}\\right)iKa_{2}^{\mu}\left(\prod_{\tilde{\lambda}=\lambda+2}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}\\right)
        && \lambda_{0}>2,\ \lambda=1,\mu=0 \\\\
        \left(\partial_{p}W_{\lambda_{0}}^{\mu}+W_{\lambda_{0}}^{\mu}\sum_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}\\frac{\partial_{p}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}{w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}}\\right)\prod_{\tilde{\lambda}=\lambda+1}^{\lambda_{0}}iKa_{\tilde{\lambda}}^{\mu}w_{\tilde{\lambda}-1}^{\mu\left(-\\right)}
        && \lambda_{0}>\lambda>1\ \mathrm{or}\ \left(\lambda_{0}>\lambda\ \mathrm{and}\ \mu\\neq0\\right)

        \\end{cases}

    """
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
        return ((dW[lam0] + W[lam0]*Sum)*(-get_a(2,0,d)/get_a(1,0,d))*Prod)
    elif lam0==2 and lam==1 and mu==0:
        return W[lam0]*dwm[1]*1j*K*get_a(2,0,d)
    elif lam0>2 and lam==1 and mu==0:
        Prod=W[lam0]*dwm[1]*1j*K*get_a(2,0,d)
        for lamt in seq(lam+2,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jm[lamt-1] 
        return Prod
    elif lam0>lam and (lam>1 or mu !=0):
        Sum=0.0
        Prod=1.0
        for lamt in seq(lam+1,lam0):
            Prod=Prod*1j*K*get_a(lamt,mu,d)/jm[lamt-1]
            Sum=Sum+dwm[lamt-1] * jm[lamt-1]
        return ((dW[lam0] + W[lam0]*Sum)*Prod) 


# --------------------------------------------------- #
# ## Residue Evaluation
# --------------------------------------------------- #



@jit(nopython=True)
def largeKResidues(K,eig,mu,nlam=10,lamMax=500,d=3):
    """
    We would like to calculate first order residues defined as

    .. math::
        Res\left(\\epsilon_{l}\\right)=\lim_{p\to\\epsilon_{l}}\left[\left(p-\\epsilon_{l}\\right)\mathcal{G}_{l}\\right]

    by L'Hospital's Rule

    .. math::
        Res\left(\\epsilon_{l}\\right)=\\frac{1}{\left.\\frac{\partial}{\partial p}\left(\\frac{1}{\mathcal{G}\left(p\\right)}\\right)\\right|_{p=\\epsilon_{l}}}

    Calculate Residues using continued fraction

    Returns a nlam x nlam numpy array of residuces with indicies (lam0,lam)

    """
    p=eig
    jp=get_jp(p,mu,K,lamMax,d)
    jm=get_jm(p,mu,K,lamMax,d)
    djp=get_djp(p,mu,K,lamMax,jp,d)
    djm=get_djm(p,mu,K,lamMax,jm,d) 
    dwp=get_dw(jp,djp,mu)
    dwm=get_dw(jm,djm,mu)
   
    invG=np.zeros((nlam,nlam),dtype=type(1+1j))
    dinvG=np.zeros((nlam,nlam),dtype=type(1+1j))
    for lam0 in range(abs(mu),nlam):
        X1=get_X1(p,lam0,mu,K,jp,jm,d)
        dX1=get_dX1(lam0,mu,K,dwp,dwm,d)
        for lam in range(abs(mu),nlam):
            if lam0==lam:
                dinvG[lam0,lam]=dX1
            elif lam0<lam:
                cp=get_cp(lam,lam0,mu,K,d)
                pm=1
                Sp=get_Sp(lam,lam0,jp,djp,pm)
            
                Xp=get_Xp(lam,lam0,jp,pm)
                dXp=get_dXp(Sp,Xp)
                
                dinvG[lam0,lam]=get_dinvG(cp,X1,dX1,Xp,dXp)
            elif lam0>lam:
                cm=get_cm(lam,lam0,mu,K,d)
                pm=-1
                Sm=get_Sp(lam,lam0,jm,djm,pm)
                
                Xm=get_Xp(lam,lam0,jm,pm)
                dXm=get_dXp(Sm,Xm)
                
                dinvG[lam0,lam]=get_dinvG(cm,X1,dX1,Xm,dXm)
    Res=1.0/dinvG # residue of 1st order pole 
    return Res



@jit(nopython=True)
def CalcInvG(K,eig,mu,nlam=10,d=3,lamMax=500):
    """Calcualte :math:`\\frac{1}{G}`

    Returns:
        list of ...

        invG (:obj:`numpy.array`): nlam x nlam numpy array with indicies(lam0,lam)

        dinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        ddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        dddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

    """
    p=eig
    jp=get_jp(p,mu,K,lamMax,d)
    jm=get_jm(p,mu,K,lamMax,d)
    djp=get_djp(p,mu,K,lamMax,jp,d)
    djm=get_djm(p,mu,K,lamMax,jm,d)
    ddjp=get_ddjp(p,mu,K,lamMax,jp,djp,d)
    ddjm=get_ddjm(p,mu,K,lamMax,jm,djm,d)
    dddjp=get_dddjp(p,mu,K,lamMax,jp,djp,ddjp,d)
    dddjm=get_dddjm(p,mu,K,lamMax,jm,djm,ddjm,d)
    
    dwp=get_dw(jp,djp,mu)
    dwm=get_dw(jm,djm,mu)
    ddwp=get_ddw(jp,djp,ddjp,mu)
    ddwm=get_ddw(jm,djm,ddjm,mu)
    dddwp=get_dddw(jp,djp,ddjp,dddjp,mu)
    dddwm=get_dddw(jm,djm,ddjm,dddjm,mu)
    
    invG=np.zeros((nlam,nlam),dtype=type(1+1j))
    dinvG=np.zeros((nlam,nlam),dtype=type(1+1j))
    ddinvG=np.zeros((nlam,nlam),dtype=type(1+1j))
    dddinvG=np.zeros((nlam,nlam),dtype=type(1+1j))
    for lam0 in range(abs(mu),nlam):
        X1=get_X1(p,lam0,mu,K,jp,jm,d)
        #dX1=X1prime(p,lam0,mu,K,jp,jm,djp,djm,d)
        dX1=get_dX1(lam0,mu,K,dwp,dwm,d)
        ddX1=get_ddX1(lam0,mu,K,ddwp,ddwm,d)
        dddX1=get_dddX1(lam0,mu,K,dddwp,dddwm,d)
        for lam in range(abs(mu),nlam):
            if lam0==lam:
                invG[lam0,lam]=X1
                dinvG[lam0,lam]=dX1
                ddinvG[lam0,lam]=ddX1
                dddinvG[lam0,lam]=dddX1
            elif lam0<lam:
                cp=get_cp(lam,lam0,mu,K,d)
                pm=1
                Sp=get_Sp(lam,lam0,jp,djp,pm)
                dSp=get_dSp(lam,lam0,jp,djp,ddjp,pm)
                ddSp=get_ddSp(lam,lam0,jp,djp,ddjp,dddjp,pm)
            
                Xp=get_Xp(lam,lam0,jp,pm)
                dXp=get_dXp(Sp,Xp)
                ddXp=get_ddXp(Sp,dSp,Xp)
                dddXp=get_dddXp(Sp,dSp,ddSp,Xp) 
                
                invG[lam0,lam]=get_invG(cp,X1,Xp)
                dinvG[lam0,lam]=get_dinvG(cp,X1,dX1,Xp,dXp)
                ddinvG[lam0,lam]=get_ddinvG(cp,X1,dX1,ddX1,Xp,dXp,ddXp)
                dddinvG[lam0,lam]=get_dddinvG(cp,X1,dX1,ddX1,dddX1,Xp,dXp,ddXp,dddXp)
            elif lam0>lam:
                cm=get_cm(lam,lam0,mu,K,d)
                pm=-1
                Sm=get_Sp(lam,lam0,jm,djm,pm)
                dSm=get_dSp(lam,lam0,jm,djm,ddjm,pm)
                ddSm=get_ddSp(lam,lam0,jm,djm,ddjm,dddjm,pm)
                
                Xm=get_Xp(lam,lam0,jm,pm)
                dXm=get_dXp(Sm,Xm)
                ddXm=get_ddXp(Sm,dSm,Xm)
                dddXm=get_dddXp(Sm,dSm,ddSm,Xm)
                
                invG[lam0,lam]=get_invG(cm,X1,Xm)
                dinvG[lam0,lam]=get_dinvG(cm,X1,dX1,Xm,dXm)
                ddinvG[lam0,lam]=get_ddinvG(cm,X1,dX1,ddX1,Xm,dXm,ddXm)
                dddinvG[lam0,lam]=get_dddinvG(cm,X1,dX1,ddX1,dddX1,Xm,dXm,ddXm,dddXm)
    return invG, dinvG, ddinvG, dddinvG



@jit(nopython=True)
def CalcG(K,p,mu,nlam=10,d=3,lamMax=500):
    """
    Calculate the propagator, :math:`G`
    """
    jp=get_jp(p,mu,K,lamMax,d)
    jm=get_jm(p,mu,K,lamMax,d)
    djp=get_djp(p,mu,K,lamMax,jp,d)
    djm=get_djm(p,mu,K,lamMax,jm,d)
    
    wp=get_w(jp,mu)
    wm=get_w(jm,mu)
    dwp=get_dw(jp,djp,mu)
    dwm=get_dw(jm,djm,mu)
    
    W=get_W(p,mu,K,lamMax,jp,jm,d=3)
    G=np.zeros((nlam,nlam),dtype=type(1.0+1j))*np.NaN
    for lam0 in range(abs(mu),nlam):
        for lam in range(abs(mu),nlam):
            G[lam0,lam]=get_G(p,lam,lam0,mu,K,jm,jp,W,d=3)
    return G, W, wp, wm




def LurentifyG(K,eig,l,mu,nlam=10,d=3,lamMax=500,
               cutoff=10**-11, 
               lowKcut=10**-7, 
               zeroKcut=10**-10,
              ): # sad face because cutoffs suck
    """
    .. math::
        R & =\lim_{p\to0}\left(p-\\epsilon_{1}\\right)G\left(p\\right) \\\\
        & =\\frac{1}{\lim_{p\to\\epsilon}\\frac{\partial}{\partial p}\left(\\frac{1}{G\left(p\\right)}\\right)}

    

    

    .. math::
        a=-\\frac{R^{2}}{2}\lim_{p\to\\epsilon}\\frac{\partial^{2}}{\partial^{2}p}\left(\\frac{1}{G}\\right)

    

    

    .. math::
        b=\\frac{a^{2}}{R}-\\frac{R^{2}}{6}\lim_{p\to\\epsilon}\\frac{\partial^{3}}{\partial^{3}p}\left(\\frac{1}{G}\\right)

    Returns:
        list of

        res (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        a (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        b (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        invG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        dinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        ddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        dddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)
    """
    if K<zeroKcut:
        return zeroKLurent(l,eig,nlam,d,mu) 
    
    old_settings = np.seterr(all='ignore')
    
    try:
        invG, dinvG, ddinvG, dddinvG = CalcInvG(K,eig,mu,nlam=nlam,d=d,lamMax=lamMax)
    except ZeroDivisionError:
        print('Encountered division by zero on CalInvG')
        print('K',K)
        print('eig',eig)
        print('mu',mu)
        raise Exception('ZeroDivisionError')
    
    tol=1.0*10**-30
    dinvG[abs(dinvG)<tol]=np.NaN
    
    res=1.0/dinvG # residue of 1st order pole
    
    # Replace with small K residue when appropriate
    for lam in range(abs(mu),nlam):
        for lam0 in range(abs(mu),nlam):
            smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
            if abs(smallK)<cutoff or K<lowKcut:
                res[lam0,lam]=smallK
    
    a=-0.5*(res**2)*ddinvG # constant term
    b=((a**2)/res) - ((res**2)/6.0)*dddinvG # (p-eps)^1 term
    # G = res/(p-eps) + a + b*(p-eps) + ...
    
    np.seterr(**old_settings)
    return res, a, b, invG, dinvG, ddinvG, dddinvG




def zeroKLurent(l,eig,nlam,d,mu):
    """
    .. math::
        \lim_{K\to0}\mathcal{G}_{\lambda\lambda_{0}}^{\mu}=\\begin{cases}

        \\frac{1}{p-\lambda\left(\lambda+d-2\\right)} &\ \ \ & \lambda=\lambda_{0} \\\\
        0 && \mathrm{otherwise}

        \\end{cases}

    

    

    .. math::
        R=\delta_{\lambda\lambda_{0}}\delta_{\lambda l}

    Returns:
        list of ...

        res (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        a (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        b (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        invG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        dinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        ddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

        dddinvG (:obj:`numpy.array`): nlam x nlam numpy array with indicies (lam0,lam)

    """
    invG=np.zeros((nlam,nlam),dtype='complex')*np.NaN
    dinvG=np.zeros((nlam,nlam),dtype='complex')
    dinvG[:,:]=np.Inf+0j
    ddinvG=np.zeros((nlam,nlam),dtype='complex')
    dddinvG=np.zeros((nlam,nlam),dtype='complex')
    a=np.zeros((nlam,nlam),dtype='complex'); a[0:mu,0:mu]=np.NaN
    b=np.zeros((nlam,nlam),dtype='complex'); b[0:mu,0:mu]=np.NaN
    res=np.zeros((nlam,nlam),dtype='complex'); res[0:mu,0:mu]=np.NaN
    
    for lam in range(mu,nlam):
        invG[lam,lam]=eig+lam*(lam+d-2)
    if abs(eig+l*(l+d-2))>10**-6:
        raise Exception('When K=0, eig sould go to l(l+d-2). l='+str(l)+' eig'+str(eig))
    
    if l<nlam:
        invG[l,l]=0.0
        dinvG[l,l]=1.0+0*1j
        res[l,l]=1.0
    
    return res, a, b, invG, dinvG, ddinvG, dddinvG




@jit(nopython=True)
def get_C(K,l,lam,mu,d=3):
    """
    .. math::
        C_{\lambda,l}^{\mu}=\\begin{cases}

        1 &\ \ \ & \lambda=l \\\\
        \left(-iK\\right)^{l-\lambda}\prod_{j=\lambda}^{l-1}\\frac{a_{j+1}^{\mu}}{l\left(l+1\\right)-j\left(j+1\\right)}
        && l>\lambda \\\\
        \left(-iK\\right)^{\lambda-l}\prod_{j=l+1}^{\lambda}\\frac{a_{j}^{\mu}}{l\left(l+1\\right)-j\left(j+1\\right)}
        && l<\lambda

        \\end{cases}

    """
    prod=1+0*1j
    Wi=l*(l+d-2)
    if l>lam:
        for j in seq(lam,l-1): # lam to l-1
            Wj=j*(j+d-2)
            prod=prod*get_a(j+1,mu,d)/(Wi-Wj)
        prod=prod*(-1j*K)**(l-lam)
    elif l<lam:
        for j in seq(l+1,lam): # l+1 to lam
            Wj=j*(j+d-2)
            prod=prod*get_a(j,mu,d)/(Wi-Wj)
        prod=prod*(-1j*K)**(lam-l)
    return prod




@jit(nopython=True)
def SmallAysmpRes(K,l,lam1,lam2,mu,d=3):
    """
    .. math::
        \lim_{p\to\\epsilon_{l}}\left[\left(p-\\epsilon_{l}\\right)\mathcal{G}_{\lambda_{0},\lambda}^{\mu}\left(K,p\\right)\\right]\\approx C_{\lambda_{0},l}^{\mu}C_{\lambda,l}^{\mu}

    l = Eigenvalue number

    lam1, lam2: spherical harmonic l index

    mu: spherical harmonic mu index

    """
    Res1=get_C(K,l,lam1,mu,d)
    Res2=get_C(K,l,lam2,mu,d)
    return Res1*Res2



@jit(nopython=True)
def residues(K,eig,l,mu,nlam=10,d=3,lamMax=500,cutoff=10**-11):
    """
    calculate residues
    using either small or large K methode as needed
    """
    # use small K limit if K<10^-3 or the limit < cutoff
    if K<10**-3:
        res=np.zeros((nlam,nlam),dtype=type(1+1j))*np.NaN
        for lam in range(abs(mu),nlam):
            for lam0 in range(abs(mu),nlam):
                #print('lam0=',lam0,'lam',lam,'mu',mu)
                smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
                res[lam0,lam]=smallK
        return res 
    
    res=largeKResidues(K,eig,mu,nlam,lamMax=lamMax)
    for lam in range(abs(mu),nlam):
        for lam0 in range(abs(mu),nlam):
            #print('lam0=',lam0,'lam',lam,'mu',mu)
            smallK=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
            if abs(smallK)<cutoff:
                res[lam0,lam]=smallK
    return res



@jit(nopython=True)
def bothResidues(K,eig,l,mu,nlam=10,d=3):
    """
    Return both large and small K residues for compairison
    """
    cutoff=10**-11
    lamMax=500
    res=largeKResidues(K,eig,mu,nlam,lamMax=lamMax)
    resSmall=res*np.NaN
    for lam in range(abs(mu),nlam):
        for lam0 in range(abs(mu),nlam):
            #print('lam0=',lam0,'lam',lam,'mu',mu)
            resSmall[lam0,lam]=SmallAysmpRes(K,l,lam,lam0,mu,d=d)
    return res, resSmall


# ============================================================ #
#
# # Autdoc
#
# ============================================================ #


import inspect
functions_to_tic= []
functions_to_auto_arg= [CalcG, CalcInvG, Epsilon, IntermediateKEigenValues, InvGprime, InvGprimeGeneral,
                         LargeKEigenValues, LurentifyG, SmallAysmpRes, X1prime, Xm_old, Xmprime, Xp_old,
                         bothResidues, get_C, get_G, get_G_zero, get_P, get_Sp, get_W, 
                         get_W_zero, get_X1, get_Xp, get_a, get_cm, get_cp, get_dG, get_dG_zero,
                         get_dSp, get_dW, get_dW_zero, get_dX1, get_dXp, get_ddSp, get_ddX1,
                         get_ddXp, get_dddX1, get_dddXp, get_dddinvG, get_dddjm, get_dddjp, get_dddw,
                         get_ddinvG, get_ddjm, get_ddjp, get_ddw, get_dinvG, get_djm, get_djm_zero,
                         get_djp, get_dw, get_dw_old, get_dwm_zero, get_invG, get_jm, get_jm_zero,
                         get_jp, get_w, get_wm_zero, largeKResidues, make_E, residues,
                         wlc_eigvals, zeroKLurent]


common_params = {'lam':'(int): :math:`l` total angular momentum enignvalue for orientation',
                 'lam1':'(int): :math:`l` total angular momentum enignvalue for orientation',
                 'lam2':'(int): :math:`l` total angular momentum enignvalue for orientation',
                 'l':'(int): eignenvalue index',
                 'mu':'(int): :math:`m` z-component of angular momentum for orientation',
                 'K':'(float): Magnitude of wavevector of fluctuation',
                 'k':'(float): Magnitude of wavevector of fluctuation',
                 'ORDEig':'(int): number of eigenvalues to include in inverse laplace transform',
                 'd':'(int): dimension e.g. 3',
                 'n_rows':'(int): number of rows.  Make this larger that ORDEig',
                 'alpha':'(float): 1/sqrt(8K)',
                 'p':'(complex): laplace inverse variable of chain length N',
                 'lamMax':'(int): maximum angular eigenvalue to include.  Higher value -> more agular detail',
                 'j':'(:obj:`numpy array`): either j(+) or j(-) as the case may be.',
                 'dj':'(:obj:`numpy array`): either dj(+) or dj(-) as the case may be.',
                 'ddj':'(:obj:`numpy array`): either ddj(+) or ddj(-) as the case may be.',
                 'dddj':'(:obj:`numpy array`): either dddj(+) or dddj(-) as the case may be.',
                 'jp':'(Iterable[complex]): j(+), Length=lamMax+1.  Obtain from get_jp ',
                 'jm':'(Iterable[complex]): j(-), Length=lamMax+1.  Obtain from get_jm ',
                 'djp':'(Iterable[complex]): dj(+)/dp, Length=lamMax+1.  Obtain from get_djp ',
                 'djm':'(Iterable[complex]): dj(-)/dp, Length=lamMax+1.  Obtain from get_djm ',
                 'ddjp':'(Iterable[complex]): dj(+)^2/d^2p, Length=lamMax+1.  Obtain from get_ddjp ',
                 'ddjm':'(Iterable[complex]): dj(-)^2/d^2p, Length=lamMax+1.  Obtain from get_ddjm ',
                 'dddjp':'(Iterable[complex]): dj(+)^3/d^3p, Length=lamMax+1.  Obtain from get_dddjp ',
                 'dddjm':'(Iterable[complex]): dj(-)^3/d^3p, Length=lamMax+1.  Obtain from get_dddjm ',
                 'wp':'(Iterable[complex]): j(+), Length=lamMax+1.  Obtain from get_w ',
                 'wm':'(Iterable[complex]): j(-), Length=lamMax+1.  Obtain from get_w ',
                 'dwp':'(Iterable[complex]): dj(+)/dp, Length=lamMax+1.  Obtain from get_dw ',
                 'dwm':'(Iterable[complex]): dj(-)/dp, Length=lamMax+1.  Obtain from get_dw ',
                 'ddwp':'(Iterable[complex]): dj(+)^2/d^2p, Length=lamMax+1.  Obtain from get_ddw ',
                 'ddwm':'(Iterable[complex]): dj(-)^2/d^2p, Length=lamMax+1.  Obtain from get_ddw ',
                 'dddwp':'(Iterable[complex]): dj(+)^3/d^3p, Length=lamMax+1.  Obtain from get_dddw ',
                 'dddwm':'(Iterable[complex]): dj(-)^3/d^3p, Length=lamMax+1.  Obtain from get_dddw ',
                 'W':'(Iterable[complex]): Infinite sum of diagrams W from equ 15 of PRE 77, 061803 (2008). Length=lamMax+1.  Obtain from get_W ',
                 'dW':'(Iterable[complex]): dW/dp, obtain from get_dW',
                 'lam0':'(int): the l spherical-harmonic component for initial orientation for the propagator.',
                 'X':'(complex): X(+) or X(-) as the case may be. From get_Xp',
                 'dX':'(complex): dX(+) or dX(-) as the case may be. From get_dXp',
                 'ddX':'(complex): ddX(+) or ddX(-) as the case may be. From get_dXp',
                 'X1':'(complex): From get_X1',
                 'Xp':'(complex): From get_Xp',
                 'Xm':'(complex): From get_Xp(...jm ...)',
                 'dX1':'(complex): From get_dX1',
                 'dXp':'(complex): From get_dXp',
                 'dXm':'(complex): From get_dXp(...jm ...)',
                 'ddX1':'(complex): From get_ddX1',
                 'ddXp':'(complex): From get_ddXp',
                 'ddXm':'(complex): From get_ddXp(...jm ...)',
                 'dddX1':'(complex): From get_dddX1',
                 'dddXp':'(complex): From get_dddXp',
                 'dddXm':'(complex): From get_dddXp(...jm ...)',
                 'cp':'(complex): Product c(+) from get_cp',
                 'cm':'(complex): Product c(-) from get_cm',
                 'S':'(complex): Sum S(+) or S(-) as the case may be.  From get_Sp',
                 'dS':'(complex): Sum dS(+) or dS(-) as the case may be.  From get_dSp',
                 'ddS':'(complex): Sum ddS(+) or ddS(-) as the case may be.  From get_ddSp',
                 'Sp':'(complex): Sum S(+) from get_Sp',
                 'Sm':'(complex): Sum S(+) from get_Sp(...jm...)',
                 'dSp':'(complex): Sum S(+) from get_dSp',
                 'dSm':'(complex): Sum S(+) from get_dSp(...jm...)',
                 'ddSp':'(complex): Sum S(+) from get_ddSp',
                 'ddSm':'(complex): Sum S(+) from get_ddSp(...jm...)',
                 'eig':'(complex): eigenvalue at which to evaluate',
                 'nlam':'(int): number of lambda values to return. nlam<lamMax',
                 'cutoff': '(float): when res<cutoff use small Aysmp',
                 'lowKcut': '(float): when to switch to low K expression for residue',
                 'zeroKcut': '(float): when to use only one residue at exactly K=0',
                 'pm': '(int): +1 or -1 as the case may be'
                }

for f in functions_to_auto_arg:
    args = list(inspect.signature(f).parameters)
    if f.__doc__ == None:
        f.__doc__='Not sure what this does \n\n    Args:'
    else:
        f.__doc__ = f.__doc__ + '\n\n    Args:'
    for arg in args:
        f.__doc__= f.__doc__ + '\n        '+arg+' '+common_params[arg]
        
f= None 

