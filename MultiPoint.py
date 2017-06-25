
# coding: utf-8

# In[ ]:

import numpy as np
import WLCgreen as wlc
import matplotlib.pyplot as plt
#%matplotlib
import propagator 
import special as sp
from numba import jit


# ## 2 point

# ###  AA

# \begin{align*}
# \sum_{l}Res_{l}\int_{0}^{N}ds_{2}\int_{0}^{s_{2}}ds_{1}e^{\epsilon_{l}\left(s_{2}-s_{1}\right)} & =\sum_{l}Res_{l}\left(\frac{e^{N\epsilon}}{\epsilon^{2}}-\frac{1}{\epsilon^{2}}-\frac{N}{\epsilon}\right)
# \end{align*}

# In[ ]:

# expl expansion methode
def IAAexplicit(N,fa,lam0,lam,p1):
    eig=p1.eig
    res=p1.res
    out=0.0
    for l in range(abs(p1.mu),p1.ORDEig):
        eps=eig[l]
        R=res[l][lam0,lam]
        out=out+R*sp.xpl(N*fa*eps,2)
    return out*((N*fa)**2) 


# \[
# \underset{AA}{I_{1}^{\left(2\right)}}\left(N\right)=\sum_{l_{1}}\frac{e^{\epsilon_{1}N}}{\epsilon_{1}^{2}}R_{1}+dG_{1}\left(0\right)+NG\left(0\right)
# \]

# In[ ]:

# residue resum methode
def IAAresum(N,fa,lam0,lam,p1):
    res=p1.res
    eig=p1.eig
    G0=p1.G0
    dG0=p1.dG0
    out = 0.0
    for l in range(abs(p1.mu),p1.ORDEig):
        out=out+np.exp(eig[l]*N*fa)*(eig[l]**-2)*res[l][lam0,lam]
        #out=out-res[l][lam0,lam]*(eig[l]**-2)
        #out=out-res[l][lam0,lam]*(N/eig[l])
        
    out=out+N*fa*G0[lam0,lam]
    out=out+dG0[lam0,lam]
    return out


# In[ ]:

# Choose best methode
def IAAswitch(N,fa,lam0,lam,p1):
    if p1.K < 10.0:
        return IAAexplicit(N,fa,lam0,lam,p1)
    else:
        return IAAresum(N,fa,lam0,lam,p1)


# ### AB

# $$
# \int_{Nf_{A}}^{N}ds_{2}\int_{0}^{Nf_{A}}ds_{1}\mathcal{G}\left(\left|s_{2}-s_{1}\right|\right)=\sum_{l=0}^{\infty}\frac{\left(e^{\epsilon_{l}Nf_{A}}-1\right)\left(e^{\epsilon_{l}N\left(1-f_{A}\right)}-1\right)}{\epsilon_{l}^{2}}R_1
# $$

# In[ ]:

# Double Intigral Methode
def IABexplicit(N,fa,lam0,lam,p1):
    eig=p1.eig
    res=p1.res
    fb=1-fa
    out=0.0
    for l in range(abs(p1.mu),p1.ORDEig):
        eps=eig[l]*N
        R=res[l][lam0,lam]
        out=out+R*sp.xpl(fa*eps,1)*sp.xpl(fb*eps,1)
    return out*(N**2)*fa*fb   


# 
# $$
# \underset{AB}{I_{1}^{\left(2\right)}}\left(N\right)=\sum_{l_{1}=0}^{\infty}R_{1}\frac{\left(e^{\epsilon_{1}N}-e^{Nf_{A}\epsilon_{1}}-e^{Nf_{B}\epsilon_{1}}\right)}{\epsilon_{1}^{2}}-dG_{1}
# $$
# 

# In[ ]:

# Resum methode
def IABresum(N,fa,lam0,lam,p1):
    eig=p1.eig
    res=p1.res
    dG0=p1.dG0
    fb=1-fa
    out=0.0
    for l in range(abs(p1.mu),p1.ORDEig):
        eps=eig[l]
        R=res[l][lam0,lam]
        out=out+R*(np.exp(N*eps) - np.exp(N*fa*eps) - np.exp(N*fb*eps))/(eps**2)
    out=out-dG0[lam0,lam]
    return out   


# In[ ]:

# Choose best methode
def IABswitch(N,fa,lam0,lam,p1):
    if p1.K < 10.0:
        return IABexplicit(N,fa,lam0,lam,p1)
    else:
        return IABresum(N,fa,lam0,lam,p1)


# ## 3 Point

# $$
# -f2resum\left(j,\epsilon_{1},\epsilon_{2},N\right)\equiv\sum_{l_{1},l_{2}}R_{1}R_{2}\frac{1}{\left(\epsilon_{1}-\epsilon_{2}\right)}\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)
# $$
# 
# 
# we may use
# 
# $$
# \sum_{l_{2}=0}^{\infty}\left(\frac{-R{}_{2}\left(\epsilon_{l}\right)}{\epsilon_{l}-p}\right)=\begin{cases}
# G_{2}\left(p\right) & p\notin\left\{ \epsilon_{2}\right\} \\
# \frac{-R{}_{2}\left(\epsilon_{2,i}\right)}{\epsilon_{2,i}-p}+a_{2,i} & p=\epsilon_{2,i}
# \end{cases}
# $$
# 
# 
# $$
# a_{2,i}=-\frac{R_{2,i}^{2}}{2}\lim_{p\to\epsilon}\frac{\partial^{2}}{\partial^{2}p}\left(\frac{1}{G_{2}\left(p\right)}\right)
# $$
# 
# 
# so
# 
# $$
# \begin{align*}
# -f2resum\left(j,N\right)= & \sum_{l_{1},\epsilon_{1}\notin\epsilon_{2}}R_{1}\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}\mathcal{G}_{2}\left(\epsilon_{1}\right)+\sum_{l_{2},\epsilon_{1}\notin\epsilon_{2}}R_{2}\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\mathcal{G}_{1}\left(\epsilon_{2}\right)\\
#  & \sum_{\epsilon_{1}=\epsilon_{2}}\left(R_{1}\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}a_{2}\left(\epsilon_{=}\right)+R_{2}\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}a_{1}\left(\epsilon_{=}\right)+\lim_{\epsilon_{1},\epsilon_{2}\to\epsilon_{=}}R_{1}R_{2}\frac{1}{\left(\epsilon_{1}-\epsilon_{2}\right)}\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)\right)
# \end{align*}
# $$

# In[ ]:

def f2resum(j,p1,p2,lam0_1,lam_1,lam0_2,lam_2,N):
    set1only, set2only, pairs      =propagator.IntersectEig2(p1.eig,p2.eig,tol=10**-6)
    out=0.0+0.0j
    
    G2at1 = p2.G_others(lam0_2,lam_2,p1)
    for l1 in range(abs(p1.mu),p1.ORDEig):
        if set1only[l1]==0:
            continue
        eps1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        out=out+(eps1**(-1*j))*np.exp(eps1*N)*R1*G2at1[l1]
        
    G1at2 = p1.G_others(lam0_1,lam_1,p2)
    for l2 in range(abs(p2.mu),p2.ORDEig):
        if set2only[l2]==0:
            continue
        eps2=p2.eig[l2]
        R2=p2.res[l2][lam0_2,lam_2]
        out=out+(eps2**(-1*j))*np.exp(eps2*N)*R2*G1at2[l2]
        
    for pair in pairs:
        (l1,l2)=pair
        eps1=p1.eig[l1]    
        eps2=p2.eig[l2] 
        R1=p1.res[l1][lam0_1,lam_1]
        R2=p2.res[l2][lam0_2,lam_2]
        a1=p1.a[l1][lam0_1,lam_1]
        a2=p2.a[l2][lam0_2,lam_2]
        out = out + R1*R2*sp.f2(j,eps1,eps2,N)
        out = out + a2*R1*np.exp(N*eps1)/(eps1**j)
        out = out + a1*R2*np.exp(N*eps2)/(eps2**j)
    
    return -out


# In[ ]:

def I3explicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,fun):
    tol = np.NaN #10**-16
    sum1=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        sum2=0.0
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e2=p2.eig[l2]
            R2=p2.res[l2][lam0_2,lam_2]
            out = R1*R2*fun(e1,e2,N,fa)
            sum2=sum2+out
            if abs(out/(sum2+tol))<tol and l2>max(lam0_2,lam_2):
                break
        sum1=sum1+sum2
        if abs(sum2/(sum1+tol))<tol and l1>max(lam0_1,lam_1):
            break
    return sum1


# ### AAA

# \begin{align*}
# I_{1,2}^{\left(3\right)}\left(N\right) & =\sum_{l_{1},l_{2}=0}^{\infty}R{}_{1}\left(\epsilon_{1}\right)R{}_{2}\left(\epsilon_{2}\right)\left(\frac{\epsilon_{1}^{-2}e^{N\epsilon_{1}}}{\epsilon_{1}-\epsilon_{2}}+\frac{\epsilon_{2}^{-2}e^{N\epsilon_{2}}}{\epsilon_{2}-\epsilon_{1}}+\frac{1}{\epsilon_{1}\epsilon_{2}^{2}}+\frac{1}{\epsilon_{1}^{2}\epsilon_{2}}+\frac{N}{\epsilon_{1}\epsilon_{2}}\right)
# \end{align*}

# In[ ]:

def IAAA(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    out=-f2resum(2,p1,p2,lam0_1,lam_1,lam0_2,lam_2,N*fa)
            
    out=out+N*fa*p2.G0[lam0_2,lam_2]*p1.G0[lam0_1,lam_1]
    out=out+p2.dG0[lam0_2,lam_2]*p1.G0[lam0_1,lam_1]
    out=out+p1.dG0[lam0_1,lam_1]*p2.G0[lam0_2,lam_2]
    return out


# \begin{align*}
# \underset{AAA}{I_{1,2}^{\left(3\right)}}\left(Nf_{A}\right)= & \sum_{l_{1},l_{2}}\left[Res_{1}\left(\epsilon_{1,l_{1}}\right)Res_{2}\left(\epsilon_{2,l_{2}}\right)f2\left(2,\epsilon_{1,l_{1}},\epsilon_{2,l_{2}},Nf_{A}\right)\right]\\
#  & +Nf_{A}\mathcal{G}_{2}\left(0\right)\mathcal{G}_{1}\left(0\right)+\mathcal{G}'_{2}\left(0\right)\mathcal{G}_{1}\left(0\right)+\mathcal{G}_{2}\left(0\right)\mathcal{G}'_{1}\left(0\right)
# \end{align*}

# In[ ]:

import warnings
def IAAAresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,debug=False):
    warnings.warn('IAAAresum is know to be inaccurate. Use IAAA')
    out=0.0+0.0j
    summ=0.0+0.0j
    for l1 in range(abs(p1.mu),p1.ORDEig):
        for l2 in range(abs(p2.mu),p2.ORDEig):
            eps1=p1.eig[l1]
            eps2=p2.eig[l2]
            R1=p1.res[l1][lam0_1,lam_1]
            R2=p2.res[l2][lam0_2,lam_2]
            
            summ=summ+sp.f2(2,eps1,eps2,N*fa)
    
    t1=N*fa*p2.G0[lam0_2,lam_2]*p1.G0[lam0_1,lam_1]
    t2=p2.dG0[lam0_2,lam_2]*p1.G0[lam0_1,lam_1]
    t3=p1.dG0[lam0_1,lam_1]*p2.G0[lam0_2,lam_2]
    out=t1+t2+t3+summ
    if debug:
        return out, t1, t2, t3, summ
    return out


# In[ ]:

def caseAAA(eps1,eps2,N,fa):
    tol=10**-6
    e1=eps1*N*fa
    e2=eps2*N*fa
    if sp.relDif(e1,e2) > tol:
        out= (sp.xpl(e2,2)-sp.xpl(e1,2))/(e2-e1)
    else:
        out= -2*sp.xpl(e1,3)+sp.xpl(e1,2)
    return out*((N*fa)**3)


# In[ ]:

def IAAAexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    return I3explicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,caseAAA)


# In[ ]:

# Choose best methode
def IAAAswitch(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    if np.sqrt(p1.K*p2.K) < 10.0/np.sqrt(N):
        return I3explicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,caseAAA)
    else:
        return IAAA(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2)


# ### ABB

# \begin{align*}
# \underset{ABB}{I_{1,2}^{\left(3\right)}}\left(Nf_{A}\right)= & \sum_{l_{1},l_{2}}Res_{1}\left(\epsilon_{1,l_{1}}\right)Res_{2}\left(\epsilon_{2,l_{2}}\right)\left[\frac{e^{Nf_{A}\epsilon_{1}}}{\epsilon_{1}}f2\left(1,\epsilon_{1},\epsilon_{2},Nf_{B}\right)-f2\left(2,\epsilon_{1,l_{1}},\epsilon_{2,l_{2}},Nf_{B}\right)\right]\\
#  & +\mathcal{G}_{1}\left(0\right)\sum_{l_{2}=0}^{\infty}Res_{2}\left(\epsilon_{2}\right)\frac{e^{Nf_{B}\epsilon_{2}}}{\epsilon_{2}^{2}}-\mathcal{G}{}_{2}\left(0\right)\sum_{l_{1}=0}^{\infty}Res_{1}\left(\epsilon_{1}\right)\frac{e^{Nf_{A}\epsilon_{1}}}{\epsilon_{1}^{2}}\\
#  & -\mathcal{G}'_{1}\left(0\right)\mathcal{G}{}_{2}\left(0\right)
# \end{align*}

# In[ ]:

def IABBresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    out=0.0+0.0j
    fb=1.0-fa
    for l1 in range(abs(p1.mu),p1.ORDEig):
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e1=p1.eig[l1]
            e2=p2.eig[l2]
            R1=p1.res[l1][lam0_1,lam_1]
            R2=p2.res[l2][lam0_2,lam_2]
            
            out=out+R1*R2*(np.exp(N*fa*e1)*sp.f2(1,e1,e2,N*fb)/e1                            - sp.f2(2,e1,e2,N*fb))
    
    temp=0.0
    for l2 in range(abs(p2.mu),p2.ORDEig):     
        e2=p2.eig[l2]
        R2=p2.res[l2][lam0_2,lam_2]
        temp=temp+R2*np.exp(N*fb*e2)/(e2**2)
    out = out + temp*p1.G0[lam0_1,lam_1]
         
    temp=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):     
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        temp=temp+R1*np.exp(N*fa*e1)/(e1**2)
    out = out - temp*p2.G0[lam0_2,lam_2]
    
    out = out - p1.dG0[lam0_1,lam_1]*p2.G0[lam0_2,lam_2]
    return out


# In[ ]:

def caseABB(eps1,eps2,N,fa):
    tol=10**-6
    e1=eps1*N
    e2=eps2*N
    fb=1.0-fa
    if sp.relDif(e1,e2) > tol:
        out= (-e1*sp.xpl(e1,3)              +fa**3*e1*sp.xpl(fa*e1,3)+fb**3*e1*sp.xpl(fb*e1,3)              +fb*fa**2*e1*sp.xpl(fa*e1,2)+fa*fb**2*e2*sp.xpl(fb*e2,2)             +(fb*fa)**2*e2*sp.xpl(fb*e2,2)*e1*sp.xpl(fa*e1,2))/(e2-e1)
    else:
        out= fa**3*sp.xpl(fa*e1,3)-fb**3*sp.xpl(fb*e1,2)             +fb**3*sp.xpl(fb*e1,3)-sp.xpl(e1,3)+fb*sp.xpl(e1,2)
    return out*(N**3)


# In[ ]:

def IABBexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    return I3explicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,caseABB)


# In[ ]:

# Choose best methode
def IABBswitch(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2):
    if np.sqrt(p1.K*p2.K) < 10.0/np.sqrt(N):
        return I3explicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2,caseABB)
    else:
        return IABBresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,p1,p2)


# ## 4 point

# ### General Functions

# $$
# \mathrm{f3resum}=\sum_{l_{1},l_{2},l_{3}}R_{1}R_{2}R_{3}\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{n}}\frac{1}{\left(\epsilon_{1}-\epsilon_{2}\right)\left(\epsilon_{3}-\epsilon_{1}\right)}+\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{n}}\frac{1}{\left(\epsilon_{1}-\epsilon_{2}\right)\left(\epsilon_{2}-\epsilon_{3}\right)}+\frac{e^{N\epsilon_{3}}}{\epsilon_{3}^{n}}\frac{1}{\left(\epsilon_{2}-\epsilon_{3}\right)\left(\epsilon_{3}-\epsilon_{1}\right)}\right)
# $$
# 
# 
# $$
# \mathrm{f3resum=\mathrm{single\ poles}+\mathrm{double\ poles}+\mathrm{triple\ poles}}
# $$
# 
# 
# $$
# \mathrm{single\ pole}=\frac{e^{\epsilon_{1}N}}{\epsilon_{1}^{2}}\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)\breve{\mathcal{G}}_{2}\left(\epsilon_{1}\right)R_{1}
# $$
# 
# 
# $$
# \mathrm{double\ pole}=\frac{e^{N\epsilon_{1}}\left(N\epsilon_{1}-2\right)}{\epsilon_{1}^{3}}R_{1}R_{2}\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)+\frac{e^{\epsilon_{1}N}}{\epsilon_{1}^{2}}\left(R_{1}R_{2}\partial\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)+a_{1}R_{2}\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)+R_{1}a_{2}\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)\right)
# $$
# 
# 
# $$
# \mathrm{triple\ pole}=\spadesuit+\diamondsuit_{23}a_{1}+\diamondsuit_{31}a_{2}+\diamondsuit_{12}a_{3}+\frac{e^{\epsilon N}}{\epsilon^{2}}\left(R_{1}a_{2}a_{3}+a_{1}R_{2}a_{3}+a_{1}a_{2}R_{3}\right)+\heartsuit
# $$
# 
# 
# $$
# \spadesuit_{123}=-R_{1}R_{2}R_{3}\frac{e^{N\epsilon}\left(N^{2}\epsilon^{2}-2nN\epsilon-n\left(n+1\right)\right)}{2\epsilon^{n+2}}
# $$
# 
# 
# $$
# \heartsuit=\frac{e^{\epsilon N}}{\epsilon^{2}}\left(2b_{1}R_{2}R_{3}+R_{1}2b_{2}R_{3}+R_{1}R_{2}2b_{3}\right)
# $$
# 
# 
# $$
# \diamondsuit_{12}=R_{1}R_{2}\frac{e^{N\epsilon}\left(N\epsilon-n\right)}{\epsilon^{n+1}}
# $$

# In[ ]:

def f3resum(j,p1,p2,p3,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3,N,fun=None):
    set1only, set2only, set3only,       pairs12, pairs23, pairs31, triplets      =propagator.IntersectEig3(p1.eig,p2.eig,p3.eig,tol=10**-10)
    
    if fun==None:
        def fun(e1):
            return 1.0
    
    out=0.0+0.0j
    
    out=out+singlePoles(set1only,p1,p2,p3,j,N,
               lam0_1,lam_1,
               lam0_2,lam_2,
               lam0_3,lam_3,fun)
        
    out=out+singlePoles(set2only,p2,p3,p1,j,N,
               lam0_2,lam_2,
               lam0_3,lam_3,
               lam0_1,lam_1,fun)
        
    out=out+singlePoles(set3only,p3,p1,p2,j,N,
               lam0_3,lam_3,
               lam0_1,lam_1,
               lam0_2,lam_2,fun)
      
    out=out+doublePoles(pairs12,p1,p2,p3,j,N,
               lam0_1,lam_1,
               lam0_2,lam_2,
               lam0_3,lam_3,fun)
        
    out=out+doublePoles(pairs23,p2,p3,p1,j,N,
               lam0_2,lam_2,
               lam0_3,lam_3,
               lam0_1,lam_1,fun)
        
    out=out+doublePoles(pairs31,p3,p1,p2,j,N,
               lam0_3,lam_3,
               lam0_1,lam_1,
               lam0_2,lam_2,fun)

    
    out=out+triplePoles(triplets,p3,p1,p2,j,N,
               lam0_3,lam_3,
               lam0_1,lam_1,
               lam0_2,lam_2,fun)
        
    out = - out # f3 is defined as -(club+club+club)
    return out


# $$
# \sum_{simple\ l_{1}}R_{1}\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{n}}\breve{\mathcal{G}}_{2}\left(\epsilon_{1}\right)\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)
# $$

# In[ ]:

def singlePoles(set1only,p1,p2,p3,j,N,
               lam0_1,lam_1,
               lam0_2,lam_2,
               lam0_3,lam_3,fun):
    G2at1 = p2.G_others(lam0_2,lam_2,p1)
    G3at1 = p3.G_others(lam0_3,lam_3,p1)
    out=0.0+0.0j
    for l1 in range(abs(p1.mu),p1.ORDEig):
        if set1only[l1]==0:
            continue
        eps1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        out=out+(eps1**(-1*j))*np.exp(eps1*N)*R1*G2at1[l1]*G3at1[l1]                *fun(eps1)  
    return out


# $$
# \sum_{double\ 12}\left[R_{1}R_{2}\lim_{\epsilon_{1}\to\epsilon_{2}}\left(\frac{\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{n}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{n}}}{\epsilon_{1}-\epsilon_{2}}\right)\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)+\frac{e^{\epsilon_{1}N}}{\epsilon_{1}^{2}}R_{1}R_{2}\partial_{p}\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)+\left(R_{2}\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{n}}a_{1}+R_{1}\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{n}}a_{2}\right)\breve{\mathcal{G}}_{3}\left(\epsilon_{1}\right)\right]
# $$

# In[ ]:

def doublePoles(pairs12,p1,p2,p3,j,N,
               lam0_1,lam_1,
               lam0_2,lam_2,
               lam0_3,lam_3,fun):
    G3at1 = p3.G_others(lam0_3,lam_3,p1)
    G3at2 = p3.G_others(lam0_3,lam_3,p2)
    dG3at1 = p3.dG_others(lam0_3,lam_3,p1)
    dG3at2 = p3.dG_others(lam0_3,lam_3,p2)
    out=0.0+0.0j
    for pair in pairs12:
        (l1,l2)=pair
        eps1=p1.eig[l1]    
        eps2=p2.eig[l2] 
        R1=p1.res[l1][lam0_1,lam_1]
        R2=p2.res[l2][lam0_2,lam_2]
        a1=p1.a[l1][lam0_1,lam_1]
        a2=p2.a[l2][lam0_2,lam_2]
        
        G3=0.5*(G3at1[l1]+G3at2[l2]) # should be about the same
        if abs(G3at1[l1]- G3at2[l2])/abs(G3) > 0.0001:
            raise Exception('G3at1[l1]- G3at2[l2] too large')
           
        eps=0.5*(eps1+eps2)
        if abs(eps1-eps2) > 0.0001:
            raise Exception('Are you sure that is a double pole?')
            
        dG3 =0.5*(dG3at1[l1]+dG3at2[l2])
        if abs(dG3at1[l1]- dG3at2[l2])/abs(dG3) > 0.0001:
            raise Exception('dG3at1[l1]+dG3at2[l2] too large') 
            
        out = out + G3*R1*R2*sp.f2(j,eps1,eps2,N)
        out = out + G3*a2*R1*np.exp(N*eps1)/(eps1**j)
        out = out + G3*a1*R2*np.exp(N*eps2)/(eps2**j)
        out = out + dG3*R1*R2*np.exp(N*eps)/(eps**j)
        
        out=out*fun(eps)
        # If we were sure eps1==eps2
        # out = out - np.exp(N*eps)*(eps**-j)*(G3*a2*R1+G3*a1*R2+dG3*R1*R2)
    return out


# In[ ]:

def triplePoles(triplets,p1,p2,p3,j,N,
               lam0_1,lam_1,
               lam0_2,lam_2,
               lam0_3,lam_3,fun):
    heart=0.0+0.0j
    Raa=0.0+0.0j
    diamond=0.0+0.0j
    spade=0.0+0.0j
    for triplet in triplets:
        (l1,l2,l3)=triplet
        eps1=p1.eig[l1]    
        eps2=p2.eig[l2] 
        eps3=p3.eig[l3]
        R1=p1.res[l1][lam0_1,lam_1]
        R2=p2.res[l2][lam0_2,lam_2]
        R3=p3.res[l3][lam0_3,lam_3]
        a1=p1.a[l1][lam0_1,lam_1]
        a2=p2.a[l2][lam0_2,lam_2]
        a3=p3.a[l3][lam0_3,lam_3]
        b1=p1.b[l1][lam0_1,lam_1]
        b2=p2.b[l2][lam0_2,lam_2]
        b3=p3.b[l3][lam0_3,lam_3] 
        
        eps = (eps1+eps2+eps3)/3.0
        
        # Raa terms
        Raa = Raa + np.exp(N*eps)/(eps**j)*(a2*a3*R1 + a3*a1*R2 + a1*a2*R3)*                    fun(eps)
        #Raa=Raa + a2*a3*R1*np.exp(N*eps1)/(eps1**j)
        #Raa=Raa + a3*a1*R2*np.exp(N*eps2)/(eps2**j)
        #Raa=Raa + a1*a2*R3*np.exp(N*eps3)/(eps3**j)
        
        # Heart term
        #heart = heart + np.exp(N*eps)/(eps**j)*(R2*R3*b1 + R3*R1*b2 + R1*R2*b3)*2
        heart=heart + R2*R3*b1*np.exp(N*eps1)/(eps1**j)*fun(eps3)
        heart=heart + R3*R1*b2*np.exp(N*eps2)/(eps2**j)*fun(eps3)
        heart=heart + R1*R2*b3*np.exp(N*eps3)/(eps3**j)*fun(eps3)
        
        # Diamond term
        diamond = diamond + np.exp(N*eps)*(N*eps-j)*(eps**(-j-1))*(R1*R2*a3 + R2*R3*a1 + R3*R1*a2)*                            fun(eps3)
        #diamond=diamond + R1*R2*sp.f2(j,eps1,eps2,N)*a3
        #diamond=diamond + R2*R3*sp.f2(j,eps2,eps3,N)*a1
        #diamond=diamond + R3*R1*sp.f2(j,eps3,eps1,N)*a2
        
        # Spade term
        spade = spade + R1*R2*R3*np.exp(N*eps)*(N**2*eps**2 - 2*j*N*eps+j*(j+1))/(2*eps**(j+2))*                        fun(eps)
        #spade=spade - R1*R2*R3*sp.f3(j,eps1,eps2,eps3,N)
        
        
    return Raa+heart+diamond+spade        


# In[ ]:

def Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,fun):
    tol = 10**-16
    sum1=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        sum2=0.0
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e2=p2.eig[l2]
            R2=p2.res[l2][lam0_2,lam_2]
            sum3=0.0
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]
                try:
                    temp = R1*R2*R3*fun(e1,e2,e3,N,fa)
                except:
                    import pdb
                    pdb.set_trace()
                sum3=sum3+temp
                if abs(temp/(sum3+tol))<tol and l3>max(lam0_3,lam_3):
                    break
            sum2=sum2+sum3
            if abs(sum3/(sum2+tol))<tol and l2>max(lam0_2,lam_2):
                break
        sum1=sum1+sum2
        if abs(sum2/(sum1+tol))<tol and l1>max(lam0_1,lam_1):
            break
    return sum1


# ### AAAA

# In[ ]:

def IAAAA(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3,K0cuttoff=10**-9):
    if p2.K < K0cuttoff:
            return IAAAAresumK2is0(N,fa,lam0_1,lam_1,                                   lam0_2,lam_2,                                   lam0_3,lam_3,p1,p3)
    
    out = - f3resum(2,p1,p2,p3,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3,N*fa)
    
    #print('f3part',out)
    
    G1=p1.G0[lam0_1,lam_1]
    G2=p2.G0[lam0_2,lam_2]
    G3=p3.G0[lam0_3,lam_3]
    dG1=p1.dG0[lam0_1,lam_1]
    dG2=p2.dG0[lam0_2,lam_2]
    dG3=p3.dG0[lam0_3,lam_3]
    zeroPole = N*fa*G1*G2*G3 + dG1*G2*G3 + G1*dG2*G3 + G1*G2*dG3
    
    #print('zeroPole',zeroPole)
    #print('ratio',out/zeroPole)
    #print('---')
    out = out + zeroPole
        
    return out


# \begin{align*}
# \underset{AAAA}{I_{1,2,3}^{\left(4\right)}\left(Nf_{A}\right)}= & -\sum_{l_{1},l_{2},l_{3}}R_{1}R_{2}R_{3}f3\left(2,\epsilon_{1},\epsilon_{2},\epsilon_{3},Nf_{A}\right)\\
#  & +Nf_{A}\mathcal{G}_{3}\mathcal{G}_{2}\mathcal{G}_{1}+\mathcal{G}'_{3}\mathcal{G}{}_{2}\mathcal{G}_{1}+\mathcal{G}_{3}\mathcal{G}'_{2}\mathcal{G}_{1}+\mathcal{G}_{3}\mathcal{G}_{2}\mathcal{G}'_{1}
# \end{align*}

# In[ ]:

@jit
def IAAAAresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3,K0cuttoff=10**-9):
    if p2.K < K0cuttoff:
            return IAAAAresumK2is0(N,fa,lam0_1,lam_1,                                   lam0_2,lam_2,                                   lam0_3,lam_3,p1,p3)
    out=0.0+0.0j
    tol = 10**-15
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e2=p2.eig[l2]
            R2=p2.res[l2][lam0_2,lam_2]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]
                #if abs(e2)<10**-7:
                #    print('epsilon of zero incountered')
                #    print('e1',e1,'e2',e2,'e3',e3)
                #    print('l1',l1,'l2',l2,'l3',l3)
                temp= - R1*R2*R3*sp.f3(2,e1,e2,e3,N*fa)
                out=out+temp
                if abs(temp/(out+tol))<tol and l3>max(lam0_3,lam_3):
                    break
    
    G1=p1.G0[lam0_1,lam_1]
    G2=p2.G0[lam0_2,lam_2]
    G3=p3.G0[lam0_3,lam_3]
    dG1=p1.dG0[lam0_1,lam_1]
    dG2=p2.dG0[lam0_2,lam_2]
    dG3=p3.dG0[lam0_3,lam_3]
    
    out = out + N*fa*G1*G2*G3 + dG1*G2*G3 + G1*dG2*G3 + G1*G2*dG3
    return out


# \[
# \sum_{l_{1},l_{3}}R{}_{1}R{}_{3}f2\left(3,\epsilon_{1},\epsilon_{3},Nf_{A}\right)+G_{1}d^{2}G_{3}+dG_{1}dG_{3}+d^{2}G_{1}G_{3}+Nf_{A}G_{1}dG_{3}+Nf_{A}dG_{1}G_{3}+\frac{N^{2}f_{A}^{2}}{2}G_{1}G_{3}
# \]

# In[ ]:

def IAAAAresumK2is0(N,fa,lam0_1,lam_1,                         lam0_2,lam_2,                         lam0_3,lam_3,p1,p3):
    if lam0_2 != lam_2:
        return 0.0
    elif lam0_2==0 and lam_2==0:
        out=0.0+0.0j
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            for l3 in range(abs(p3.mu),p3.ORDEig): 
                e3=p1.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]

                out=out+R1*R3*sp.f2(3,e1,e3,N*fa)

        G1=p1.G0[lam0_1,lam_1]
        G3=p3.G0[lam0_3,lam_3]
        dG1=p1.dG0[lam0_1,lam_1]
        dG3=p3.dG0[lam0_3,lam_3]    

        #This look can be replaced by next derivative
        ddG1=0.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            ddG1=ddG1-R1/(e1**3) 

        #This look can be replaced by next derivative
        ddG3=0.0
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]
            ddG3=ddG3-R3/(e3**3) 

        out=out+G1*ddG3 + dG1*dG3 + ddG1*G3 +             N*fa*G1*dG3 + N*fa*dG1*G3+0.5*N*fa*N*fa*G1*G3
        return out
    else: # lam_2 == lam_2
        out=0.0+0.0j
        tol = 10**-15
        e2=-lam_2*(lam_2 + p1.d -2)
        R2=1.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]
                temp= - R1*R2*R3*sp.f3(2,e1,e2,e3,N*fa)
                out=out+temp

        
        G1=p1.G0[lam0_1,lam_1]
        G2=-1.0/e2
        G3=p3.G0[lam0_3,lam_3]
        dG1=p1.dG0[lam0_1,lam_1]
        dG2=-1.0/(e2**2)
        dG3=p3.dG0[lam0_3,lam_3]

        out = out + N*fa*G1*G2*G3 + dG1*G2*G3 + G1*dG2*G3 + G1*G2*dG3        
        return out
    


# For IAAAA if $\epsilon_{1}\not\approx\epsilon_{2}\not\approx\epsilon_{3}$
# 
# \[
# =\frac{\left(\frac{\epsilon_{2}-\epsilon_{3}}{\epsilon_{1}^{2}}\right)\mathrm{expl}\left(4,N\epsilon_{1}\right)+\left(\frac{\epsilon_{3}-\epsilon_{1}}{\epsilon_{2}^{2}}\right)\mathrm{expl}\left(4,N\epsilon_{2}\right)+\left(\frac{\epsilon_{1}-\epsilon_{2}}{\epsilon_{3}^{2}}\right)\mathrm{expl}\left(4,N\epsilon_{3}\right)}{\left(\epsilon_{1}-\epsilon_{2}\right)\left(\epsilon_{1}-\epsilon_{3}\right)\left(\epsilon_{2}-\epsilon_{3}\right)}
# \]
# 
# 
# if $\epsilon_{1}\not\approx\epsilon_{2}\approx\epsilon_{3}$ 
# 
# \[
# =\frac{-N\epsilon_{2}\epsilon_{1}^{3}\mathrm{expl\left(3,N\epsilon_{2}\right)+2\epsilon_{1}^{3}\mathrm{expl}\left(4,N\epsilon_{2}\right)+\epsilon_{2}^{3}\mathrm{expl}}\left(4,N\epsilon_{1}\right)-3\epsilon_{1}^{2}\epsilon_{2}\mathrm{expl}\left(4,N\epsilon_{2}\right)+N\epsilon_{2}^{2}\epsilon_{1}^{2}\mathrm{expl}\left(3,N\epsilon_{2}\right)}{\epsilon_{1}^{4}\epsilon_{2}^{3}-2\epsilon_{1}^{3}\epsilon_{2}^{4}+\epsilon_{1}^{2}\epsilon_{2}^{5}}
# \]
# 
# 
# if $\epsilon_{1}\approx\epsilon_{2}\approx\epsilon_{3}$
# 
# \[
# =\frac{-6\mathrm{expl}\left(4,N\epsilon_{1}\right)-N^{2}\epsilon_{1}^{2}\mathrm{expl}\left(2,N\epsilon_{1}\right)+4N\epsilon_{1}\mathrm{expl}\left(3,N\epsilon_{1}\right)}{2\epsilon_{1}^{4}}
# \]

# In[ ]:

@jit(nopython=True)
def caseAAAA(eps1,eps2,eps3,N,fa):
    e3,e2,e1 = sp.arrange(eps1*N*fa,eps2*N*fa,eps3*N*fa)
    tol=10**-6
    if sp.relDif(e2,e3)>tol:
        out=-( (e2-e3)*e1**2*sp.xpl(e1,4) +                (e3-e1)*e2**2*sp.xpl(e2,4) +                (e1-e2)*e3**2*sp.xpl(e3,4) )/              ((e1-e2)*(e2-e3)*(e3-e1))
    elif min(sp.relDif(e1,e2),sp.relDif(e3,e1)) > tol:
        out=(-e2*e1*sp.xpl(e2,3)+2*e2*e1*sp.xpl(e2,4)+             e1**2*sp.xpl(e1,4)-3*e2**2*sp.xpl(e2,4)+e2**2*sp.xpl(e2,3))             /((e1-e2)**2)
    else:
        out = 3*sp.xpl(e1,4)+0.5*sp.xpl(e1,2)-2*sp.xpl(e1,3)
    return out*((N*fa)**4)


# In[ ]:

def IAAAAexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3):
    return Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAAAA)


# In[ ]:

# Use explicit calculation for low K and IAAAA for high K
def IAAAAswitch(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3):
    K=(p1.K*p2.K*p3.K)**(1/3.0)
    if K<50.0:
        out = Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAAAA)
    else:
        out =IAAAA(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3)
    return out


# ### AAAB

# In[ ]:

def IAAAB(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3,K0cuttoff=10**-9):
    raise ValueError('IAAAB not working. Use IAAABresum.')
    if p2.K < K0cuttoff:
        return IAAABresumK2is0(N,fa,lam0_1,lam_1                               ,lam0_2,lam_2                               ,lam0_3,lam_3,p1,p3)
    fb=1.0-fa
    out=0.0+0.0j
    term1=0.0+0.0j
    tol = 10**-15
    def fun(eps):
        return np.exp(eps*N*fb)/eps
    
    out = f3resum(2,p1,p2,p3,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3,N)
    out = out - f3resum(1,p1,p2,p3,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3,N,fun)

    G1=p1.G0[lam0_1,lam_1]
    G2=p2.G0[lam0_2,lam_2]
    G3=p3.G0[lam0_3,lam_3]          
                
    temp=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e1=p1.eig[l1]
            e2=p2.eig[l2]
            R1=p1.res[l1][lam0_1,lam_1]
            R2=p2.res[l2][lam0_2,lam_2]   
            temp=temp+R1*R2*sp.f2(2,e1,e2,N*fa)
    out=out+temp*G3
            
    temp=0.0    
    for l3 in range(abs(p3.mu),p3.ORDEig):
        e3=p3.eig[l3]
        R3=p3.res[l3][lam0_3,lam_3]
        temp=temp+R3*np.exp(N*fb*e3)/(e3**2)
    out=out-temp*G1*G2
    
    out = out - G1*G2*p3.dG0[lam0_3,lam_3] 
    return out


# \begin{align*}
# \underset{AAAB}{I_{1,2,3}^{\left(4\right)}\left(N\right)}= & \sum_{l_{1},l_{2},l_{3}}R{}_{1}R{}_{2}R{}_{3}\left[f3\left(2,\epsilon_{1},\epsilon_{2},\epsilon_{3},Nf_{A}\right)-\frac{e^{Nf_{B}\epsilon_{3}}}{\epsilon_{3}}f3\left(1,\epsilon_{1},\epsilon_{2},\epsilon_{3},Nf_{A}\right)\right]\\
#  & +\mathcal{G}_{3}\sum_{l_{1},l_{2}}R{}_{1}R{}_{2}f2\left(2,\epsilon_{1},\epsilon_{2},Nf_{A}\right)\\
#  & -\mathcal{G}_{1}\mathcal{G}_{2}\sum_{l_{3}}R_{3}\frac{e^{Nf_{B}\epsilon_{3}}}{\epsilon_{3}^{2}}\\
#  & -\mathcal{G}_{1}\mathcal{G}_{2}\mathcal{G}'_{3}
# \end{align*}

# In[ ]:

#@jit#(nopython=True)
def IAAABresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3,K0cuttoff=10**-9):
    if p2.K < K0cuttoff:
        return IAAABresumK2is0(N,fa,lam0_1,lam_1                               ,lam0_2,lam_2                               ,lam0_3,lam_3,p1,p3)
    fb=1.0-fa
    out=0.0+0.0j
    tol = 10**-15
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e2=p2.eig[l2]
            R2=p2.res[l2][lam0_2,lam_2]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]

                temp= R1*R2*R3*( sp.f3(2,e1,e2,e3,N*fa) -                            np.exp(N*fb*e3)*sp.f3(1,e1,e2,e3,N*fa)/e3 )
                out=out+temp
                if abs(temp/(out+tol))<tol and l3>max(lam0_3,lam_3):
                    break

    G1=p1.G0[lam0_1,lam_1]
    G2=p2.G0[lam0_2,lam_2]
    G3=p3.G0[lam0_3,lam_3]          
                
    temp=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e1=p1.eig[l1]
            e2=p2.eig[l2]
            R1=p1.res[l1][lam0_1,lam_1]
            R2=p2.res[l2][lam0_2,lam_2]   
            temp=temp+R1*R2*sp.f2(2,e1,e2,N*fa)
    out=out+temp*G3
            
    temp=0.0    
    for l3 in range(abs(p3.mu),p3.ORDEig):
        e3=p3.eig[l3]
        R3=p3.res[l3][lam0_3,lam_3]
        temp=temp+R3*np.exp(N*fb*e3)/(e3**2)
    out=out-temp*G1*G2
    
    out = out - G1*G2*p3.dG0[lam0_3,lam_3] 
    return out


# \begin{align*}
# \underset{AAAB}{I_{1,2,3}^{\left(4\right)}\left(N\right)}= & \sum_{l_{1},l_{3}}R{}_{1}R{}_{3}f2\left(2,\epsilon_{1}\epsilon_{3},Nf_{A}\right)\frac{\exp\left(Nf_{B}\epsilon_{3}-1\right)}{\epsilon_{3}}\\
#  & -G_{1}\sum_{l_{3}}R{}_{3}\frac{e^{Nf_{B}\epsilon_{3}}\left(1+Nf_{A}\epsilon_{3}\right)}{\epsilon_{3}^{3}}-dG_{1}\sum_{l_{3}}R{}_{3}\frac{e^{Nf_{B}\epsilon_{3}}}{\epsilon_{3}^{2}}\\
#  & -G_{1}d^{2}G_{3}-dG_{1}dG_{3}-Nf_{A}G_{1}dG_{3}
# \end{align*}

# In[ ]:

def IAAABresumK2is0(N,fa,lam0_1,lam_1                        ,lam0_2,lam_2                        ,lam0_3,lam_3                        ,p1,p3):
    if lam0_2 != lam_2:
        return 0.0
    elif lam0_2==0 and lam_2==0:
        fb=1.0-fa
        out=0.0+0.0j
        for l1 in range(abs(p1.mu),p1.ORDEig):
            for l3 in range(abs(p3.mu),p3.ORDEig): 
                e1=p1.eig[l1]
                e3=p1.eig[l3]
                R1=p1.res[l1][lam0_1,lam_1]
                R3=p3.res[l3][lam0_3,lam_3]

                out=out+R1*R3*sp.f2(2,e1,e3,N*fa)*(np.exp(N*fb*e3)-1)/e3

        G1=p1.G0[lam0_1,lam_1]
        dG1=p1.dG0[lam0_1,lam_1]
        dG3=p3.dG0[lam0_3,lam_3]   

        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]
            out=out-G1*R3*np.exp(N*fb*e3)*(1+N*fa*e3)/(e3**3)
            out=out-dG1*R3*np.exp(N*fb*e3)/(e3**2)

        #This look can be replaced by next derivative
        ddG3=0.0
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]
            ddG3=ddG3-R3/(e3**3) 

        out=out-G1*ddG3-dG1*dG3-N*fa*G1*dG3
        return out
    else: # lam_2 == lam_2
        fb=1.0-fa
        out=0.0+0.0j
        tol = 10**-15
        e2=-lam_2*(lam_2 + p1.d -2)
        R2=1.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]

                temp= R1*R2*R3*( sp.f3(2,e1,e2,e3,N*fa) -                                 np.exp(N*fb*e3)*                                 sp.f3(1,e1,e2,e3,N*fa)/e3 )
                out=out+temp

        G1=p1.G0[lam0_1,lam_1]
        G2=-1.0/e2
        G3=p3.G0[lam0_3,lam_3]          

        temp=0.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]  
            temp=temp+R1*R2*sp.f2(2,e1,e2,N*fa)
        out=out+temp*G3

        temp=0.0    
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]
            temp=temp+R3*np.exp(N*fb*e3)/(e3**2)
        out=out-temp*G1*G2

        out = out - G1*G2*p3.dG0[lam0_3,lam_3] 
        
        return out        


# In[ ]:

@jit(nopython=True)
def AAABwhen3equal2not1(e1,e2,fa):
    fb=1.0-fa
    return (fa*(sp.xpl(fa*e1,4)*e1**3 - sp.xpl(fa*e2,4)*e2**3)/         ((e1-e2)**2)         +fa*(e2**2*sp.xpl(fa*e2,4)-e2**2*sp.xpl(fa*e2,3))/(e1-e2)         +1.0/6.0)


# In[ ]:

@jit(nopython=True)
def caseAAAB(eps1,eps2,eps3,N,fa):
    e1,e2,e3=eps1*N,eps2*N,eps3*N
    tol=10**-6
    fb=1-fa
    ab=sp.relDif(e1,e2)<tol
    bc=sp.relDif(e2,e3)<tol
    ca=sp.relDif(e3,e1)<tol
    if (not ab) and (not bc) and (not ca):
        out=fb*fa**3*sp.xpl(fb*e3,1)*        (fa*((e2-e3)*sp.xpl(fa*e1,4)*e1**3            +(e3-e1)*sp.xpl(fa*e2,4)*e2**3            +(e1-e2)*sp.xpl(fa*e3,4)*e3**3 )/         ((e1-e2)*(e1-e3)*(e2-e3))+1.0/6.0)
    elif (not ab) and (not ca): # e2=e3
        out = fb*fa**3*sp.xpl(fb*e3,1)*AAABwhen3equal2not1(e1,e2,fa)
    elif (not ab) and (not bc): # e1=e3
        out = fb*fa**3*sp.xpl(fb*e3,1)*AAABwhen3equal2not1(e2,e3,fa)
    elif (not bc) and (not ca): # e1=e2
        out = fb*fa**3*sp.xpl(fb*e3,1)*AAABwhen3equal2not1(e3,e1,fa)
    else:
        out = sp.xpl(e1,4)-fb**4*sp.xpl(fb*e1,4)-fa**4*sp.xpl(fa*e1,4)              + 0.5*fa**2*sp.xpl(e1,2)-0.5*fa**4*sp.xpl(fa*e1,2)              +fa**4*sp.xpl(fa*e1,3)-fa*sp.xpl(e1,3)
    return out*(N**4)


# In[ ]:

def IAAABexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3):
    return Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAAAB)


# In[ ]:

def IAAABswitch(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3):
    K=(p1.K*p2.K*p3.K)**(1.0/3.0)
    if K<10:
        out = Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAAAB)
    else:
        out =IAAABresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3)
    return out


# ### AABB

# \begin{align*}
# \underset{AABB}{I_{1,2,3}^{\left(4\right)}\left(N\right)}= & \sum_{l_{1},l_{2},l_{3}}R{}_{1}R{}_{2}R{}_{3}f\left(1,\epsilon_{1},\epsilon_{2},Nf_{A}\right)f\left(1,\epsilon_{2},\epsilon_{3},Nf_{B}\right)\\
#  & -\mathcal{G}_{3}\sum_{l_{1},l_{2}}R{}_{1}R{}_{2}f\left(2,\epsilon_{1},\epsilon_{2},Nf_{A}\right)-\mathcal{G}_{1}\sum_{l_{2},l_{3}}R{}_{2}R{}_{3}f\left(2,\epsilon_{2},\epsilon_{3},Nf_{B}\right)\\
#  & +\mathcal{G}_{2}\mathcal{G}_{3}\sum_{l_{1}}R_{1}\frac{e^{Nf_{A}\epsilon_{1}}}{\epsilon_{1}^{2}}+\mathcal{G}_{1}\mathcal{G}_{2}\sum_{l_{3}}R_{3}\frac{e^{Nf_{B}\epsilon_{3}}}{\epsilon_{3}^{2}}\\
#  & -\mathcal{G}_{1}\mathcal{G}'_{2}\mathcal{G}_{3}
# \end{align*}

# In[ ]:

def IAABBresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3,K0cuttoff=10**-9):
    if p2.K < K0cuttoff:
            return IAABBresumK2is0(N,fa,lam0_1,lam_1                                       ,lam0_2,lam_2                                       ,lam0_3,lam_3,p1,p3)
    fb=1.0-fa
    out=0.0+0.0j
    tol = 10**-15
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e2=p2.eig[l2]
            R2=p2.res[l2][lam0_2,lam_2]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]

                temp= R1*R2*R3*( sp.f2(1,e1,e2,N*fa)*sp.f2(1,e2,e3,N*fb) )
                out=out+temp
                if abs(temp/(out+tol))<tol and l3>max(lam0_3,lam_3):
                    break
    G1=p1.G0[lam0_1,lam_1]
    G2=p2.G0[lam0_2,lam_2]
    G3=p3.G0[lam0_3,lam_3]
           
    temp=0.0
    for l1 in range(abs(p1.mu),p1.ORDEig):
        for l2 in range(abs(p2.mu),p2.ORDEig):
            e1=p1.eig[l1]
            e2=p2.eig[l2]
            R1=p1.res[l1][lam0_1,lam_1]
            R2=p2.res[l2][lam0_2,lam_2]   
            temp=temp+R1*R2*sp.f2(2,e1,e2,N*fa)
    out=out-temp*G3

    temp=0.0
    for l2 in range(abs(p2.mu),p2.ORDEig):
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            e2=p2.eig[l2]
            R3=p3.res[l3][lam0_3,lam_3]
            R2=p2.res[l2][lam0_2,lam_2]   
            temp=temp+R2*R3*sp.f2(2,e2,e3,N*fb)
    out=out-temp*G1    

    temp=0.0    
    for l1 in range(abs(p1.mu),p1.ORDEig):
        e1=p1.eig[l1]
        R1=p1.res[l1][lam0_1,lam_1]
        temp=temp+R1*np.exp(N*fa*e1)/(e1**2)
    out=out+temp*G2*G3    
    
    temp=0.0    
    for l3 in range(abs(p3.mu),p3.ORDEig):
        e3=p3.eig[l3]
        R3=p3.res[l3][lam0_3,lam_3]
        temp=temp+R3*np.exp(N*fb*e3)/(e3**2)
    out=out+temp*G1*G2
    
    out = out - G1*p2.dG0[lam0_2,lam_2]*G3
    return out


# \[
# \underset{AABB}{I_{1,2,3}^{\left(4\right)}\left(N\right)}=\underset{AA}{I_{1}^{\left(2\right)}}\left(N\right)\underset{BB}{I_{3}^{\left(2\right)}}\left(N\right)
# \]

# In[ ]:

def IAABBresumK2is0(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                    ,p1,p3):
    if lam0_2 != lam_2:
        return 0.0
    elif lam0_2==0 and lam_2==0:
        fb=1.0-fa
        out=IAAresum(N,fa,lam0_1,lam_1,p1)*IAAresum(N,fb,lam0_3,lam_3,p3)
        return out
    else:
        fb=1.0-fa
        out=0.0+0.0j
        tol = 10**-15
        e2=-lam_2*(lam_2 + p1.d -2)
        R2=1.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            for l3 in range(abs(p3.mu),p3.ORDEig):
                e3=p3.eig[l3]
                R3=p3.res[l3][lam0_3,lam_3]

                temp= R1*R2*R3*( sp.f2(1,e1,e2,N*fa)*sp.f2(1,e2,e3,N*fb) )
                out=out+temp
                
        G1=p1.G0[lam0_1,lam_1]
        G2=-1.0/e2
        G3=p3.G0[lam0_3,lam_3]

        temp=0.0
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]  
            temp=temp+R1*R2*sp.f2(2,e1,e2,N*fa)
        out=out-temp*G3

        temp=0.0
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]  
            temp=temp+R2*R3*sp.f2(2,e2,e3,N*fb)
        out=out-temp*G1    

        temp=0.0    
        for l1 in range(abs(p1.mu),p1.ORDEig):
            e1=p1.eig[l1]
            R1=p1.res[l1][lam0_1,lam_1]
            temp=temp+R1*np.exp(N*fa*e1)/(e1**2)
        out=out+temp*G2*G3    

        temp=0.0    
        for l3 in range(abs(p3.mu),p3.ORDEig):
            e3=p3.eig[l3]
            R3=p3.res[l3][lam0_3,lam_3]
            temp=temp+R3*np.exp(N*fb*e3)/(e3**2)
        out=out+temp*G1*G2

        out = out - G1*(-1.0/(e2**2))*G3
        return out


# In[ ]:

@jit(nopython=True)
def cow(e1,e2,fa):
    return (sp.xpl(fa*e1,2)*e1 - sp.xpl(fa*e2,2)*e2)/(e1-e2)
@jit(nopython=True)
def bull(e1,fa):
    return sp.xpl(fa*e1,1) - sp.xpl(fa*e1,2)


# In[ ]:

@jit(nopython=True)
def caseAABB(eps1,eps2,eps3,N,fa):
    tol=10**-6
    fb=1-fa
    e1,e2,e3=eps1*N,eps2*N,eps3*N
    ab=sp.relDif(e1,e2)<tol
    bc=sp.relDif(e2,e3)<tol
    ca=sp.relDif(e3,e1)<tol
    if (not ab) and (not bc) and (not ca):
        out = cow(e1,e2,fa)*cow(e2,e3,fb)
    elif (not ab) and (not ca): # e2=e3
        out = cow(e1,e2,fa)*bull(e2,fb)
    elif (not ab) and (not bc): # e1=e3
        out = cow(e1,e2,fa)*cow(e2,e3,fb)
    elif (not bc) and (not ca): # e1=e2
        out = bull(e1,fa)*cow(e2,e3,fb)
    else:
        out = bull(e1,fa)*bull(e1,fb)
    return out*(N**4)*fa*fa*fb*fb


# In[ ]:

def IAABBexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3):
    return Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAABB)


# In[ ]:

def IAABBswitch(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3):
    K=(p1.K*p2.K*p3.K)**(1/3.0)
    if K<10:
        out = Iexplicit(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3                      ,p1,p2,p3,caseAABB)
    else:
        out =IAABBresum(N,fa,lam0_1,lam_1,lam0_2,lam_2,lam0_3,lam_3               ,p1,p2,p3)
    return out

