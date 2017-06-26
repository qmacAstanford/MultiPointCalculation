
# coding: utf-8

# Here are some functions that are useful and have some numerical limits that have been fixed.

# In[ ]:

import numpy as np
from numba import jit


# \[
# \mathrm{expl\left(x,n\right)=e^{x}-\sum_{j=0}^{n-1}\frac{x^{j}}{j!}}
# \]
# 
# 
# \[
# \mathrm{expl\left(x,n\right)=\sum_{j=n}^{\infty}\frac{x^{j}}{j!}}
# \]

# In[ ]:

#@jit(nopython=True)
def factorial():
    out=int(1)
    fact = np.ones(100)
    for ii in range(1,100):
        out=out*ii
        fact[ii] = out
    return fact
fact=factorial()


# In[ ]:

@jit(nopython=True)
def expl(x,n,maxn=25):
    x=x*1.0
    if abs(x)<1.0:
        out=0.0
        for ii in range(n,maxn):
            out=out+x**ii/float(fact[ii])#float(np.math.factorial(ii))
    else:
        out=np.exp(x)
        for ii in range(0,n):
            out=out-x**ii/float(fact[ii])#float(np.math.factorial(ii))
    return out


# In[ ]:

@jit(nopython=True)
def xpl(x,n,maxn=25):
    x=x*1.0
    if abs(x)<1.0:
        out=0.0
        for ii in range(n,maxn):
            out=out+x**(ii-n)/float(fact[ii])#float(np.math.factorial(ii))
    else:
        out=np.exp(x)/(x**n)
        for ii in range(0,n):
            out=out-x**(ii-n)/float(fact[ii])#float(np.math.factorial(ii))
    return out


# In[ ]:

# sequence of integers from first to last
# for example seq(2,5)=2,3,4,5
@jit(nopython=True)
def seq(first,last):
    if first<=last:
        return np.arange(first,last+1)
    else:
        return np.arange(first,last-1,-1)


# \[
# f2\left(j,\epsilon_{1},\epsilon_{2},N\right)\equiv\frac{1}{\left(\epsilon_{1}-\epsilon_{2}\right)}\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)
# \]
# 
# 
# This is not numerically well behaved if $\epsilon_{1}\approx\epsilon_{2}$
# in which case we use:
# 
# \[
# x\equiv\frac{\epsilon_{2}-\epsilon_{1}}{2},\ \ y\equiv\frac{\epsilon_{2}+\epsilon_{1}}{2}
# \]
# 
# 
# \[
# f2=\begin{cases}
# \sum_{n=\mathrm{even}}^{\infty}x^{n}\frac{e^{Ny}}{y^{n+2}}\sum_{m=0}^{n+1}-\frac{\left(-Ny\right)^{m}}{m!} & \mathrm{if}\ j=1\\
# \sum_{n=\mathrm{even}}^{\infty}x^{n}\frac{e^{Ny}}{y^{n+3}}\sum_{m=0}^{n+1}-\frac{\left(-Ny\right)^{m}\left(n+2-m\right)}{m!} & \mathrm{if}\ j=2\\
# \sum_{n=\mathrm{even}}^{\infty}x^{n}\frac{e^{Ny}}{y^{n+4}}\sum_{m=0}^{n+1}-\frac{\left(-Ny\right)^{m}\left(n+3-m\right)\left(n+2-m\right)}{m!} & \mathrm{if}\ j=3
# \end{cases}
# \]
# 
# 
# if we perform the sum to $\infty$ this is in principle exact for
# all imput.
# 
# In order to determine whether $\epsilon_{1}\approx\epsilon_{2}$ a
# good method cutoff to use is
# 
# \[
# \frac{\left|\epsilon_{1}-\epsilon_{2}\right|}{\mathrm{min}\left(\left|\epsilon_{1}\right|,\left|\epsilon_{2}\right|\right)}=0.001
# \]
# 

# In[ ]:

@jit(nopython=True)
def f2(j,e1,e2,N):
    tol=0.001
    nmax=15
    tol=10**-15 # How may digits of accuracy 
    if relDif(e1,e2)<tol:
        x=(e2-e1)/2.0
        y=(e2+e1)/2.0
        if j == 1:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-(-N*y)**m / float(fact[m])
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-2))
                if n>4 and abs(temp) < abs(out)*tol:
                    return out
        elif j == 2:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-(n+2-m)*(-N*y)**m / float(fact[m])
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-3))
                if n>4 and abs(temp) < abs(out)*tol:
                    return out
        elif j== 3:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-0.5*(n+3-m)*(n+2-m)*(-N*y)**m / float(fact[m])
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-4))
                if n>4 and abs(temp) < abs(out)*tol:
                    return out
    else:
        out=(np.exp(N*e1)/(e1**j) - np.exp(N*e2)/(e2**j))/(e1-e2)
    return out


# In[ ]:

# Orders a, b, and c such that c is furtest from the other two
@jit(nopython=True)
def arrange(a,b,c):
    ab = relDif(a,b)
    bc = relDif(b,c)
    ac = relDif(a,c)
    if ab <= min(bc,ac):
        return (a,b,c)
    elif ac <= min(ab,bc):
        return (a,c,b)
    else:
        return (b,c,a)


# In[ ]:

@jit(nopython=True)
def relDif(a,b):
    smaller=min(abs(a),abs(b))
    if smaller <= 10**-200:
        if abs(a-b)<= 10**-200:
            return 0.0
        else:
            return np.inf
    return abs(a-b)/smaller


# \[
# f3\left(j,\epsilon_{1},\epsilon_{2},\epsilon_{3}\right)\equiv\frac{\left(\frac{e^{N\epsilon_{3}}}{\epsilon_{3}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)\epsilon_{1}+\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{3}}}{\epsilon_{3}^{j}}\right)\epsilon_{2}+\left(\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}-\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}\right)\epsilon_{3}}{\left(\epsilon_{1}-\epsilon_{2}\right)\left(\epsilon_{2}-\epsilon_{3}\right)\left(\epsilon_{3}-\epsilon_{1}\right)}
# \]
# 
# 
# in the even that $\epsilon_{1}\approx\epsilon_{2}\neq\epsilon_{3}$
# we use the alternative calculation method:
# 
# \[
# f3=\frac{\epsilon_{3}^{-j}e^{N\epsilon_{3}}+\epsilon_{2}\epsilon_{1}f2\left(j+1,\epsilon_{1},\epsilon_{2},N\right)-\epsilon_{3}f2\left(j,\epsilon_{1},\epsilon_{2},N\right)}{\left(\epsilon_{2}-\epsilon_{3}\right)\left(\epsilon_{3}-\epsilon_{1}\right)}
# \]
# 
# 
# and when $\epsilon_{1}\approx\epsilon_{2}\approx\epsilon_{3}$ we
# use the following calculation method:
# 
# \[
# f3=\begin{cases}
# \sum_{k=0}^{\infty}\frac{e^{N\epsilon_{1}}}{\left(-\epsilon_{1}\right)^{k+3}}\left(\sum_{m=0}^{k+2}\frac{\left(-N\epsilon_{1}\right)^{m}}{m!}\right)\left(\sum_{m=0}^{k}\left(\epsilon_{2}-\epsilon_{1}\right)^{k-m}\left(\epsilon_{3}-\epsilon_{1}\right)^{m}\right) & \mathrm{if\ }j=1\\
# \sum_{k=0}^{\infty}\frac{\left(-1\right)^{k+1}e^{N\epsilon_{1}}}{\epsilon_{1}^{k+4}}\left(\sum_{m=0}^{k+2}\left(-N\epsilon_{1}\right)^{m}\frac{\left(k-m+3\right)}{m!}\right)\left(\sum_{m=0}^{k}\left(\epsilon_{2}-\epsilon_{1}\right)^{k-m}\left(\epsilon_{3}-\epsilon_{1}\right)^{m}\right) & \mathrm{if\ }j=2
# \end{cases}
# \]
# 
# 
# This expression works best when $\epsilon_{1}$ is chosen to be the
# $\epsilon$ value furthest from zero. If the sum fails to converge
# revert to the previous method of calculation.
# \end{document}

# In[ ]:

@jit(nopython=True)
def f3(j,eps1,eps2,eps3,N):
    tol=0.001
    e1,e2,e3 = arrange(eps1,eps2,eps3)
    if relDif(e1,e2) >= tol: # none close to eachother
        return ((np.exp(N*e3)/(e3**j) - np.exp(N*e2)/(e2**j))*e1+                 (np.exp(N*e1)/(e1**j) - np.exp(N*e3)/(e3**j))*e2+                 (np.exp(N*e2)/(e2**j) - np.exp(N*e1)/(e1**j))*e3)/                 ((e1-e2)*(e2-e3)*(e3-e1))
    elif( max(relDif(e1,e2),relDif(e2,e3),relDif(e3,e1))  <= tol):
        # all close to eachother
        
        # order so that e1 is furthest from zero  
        # and e3 is closest to zero
        if abs(e1)<1.0:
            if abs(e2) >= abs(e1) and abs(e2)>= abs(e3):
                e1, e2, e3 = e2, e1, e3
            elif abs(e3) >= abs(e1) and abs(e3) >= abs(e2):
                e1, e2, e3 = e3, e2, e1
            if abs(e2)<abs(e3):
                e1, e2, e3 = e1, e3, e2
        kmax=20 # How many terms to try for convergence
        tol=10**-12 # How may digits of accuracy 
        if j==1:
            out=0.0
            for k in range(0,kmax):
                prod=np.exp(N*e1)*(-e1)**(-k-3)
                summ=0.0
                for m in range(0,k+2+1):
                    summ=summ+(-N*e1)**m/fact[m]
                prod=prod*summ
                summ=0.0
                for m in range(0,k+1):
                    summ=summ+(e2-e1)**(k-m) * (e3-e1)**m
                temp=prod*summ
                out=out + temp
                if k>4 and abs(temp) < abs(out)*tol:
                    return out
        elif j==2:
            out=0.0
            for k in range(0,kmax):
                prod=-np.exp(N*e1)*(-e1)**(-k-4)
                summ=0.0
                for m in range(0,k+2+1):
                    summ=summ+(-N*e1)**m *(k-m+3)/fact[m]
                prod=prod*summ
                summ=0.0
                for m in range(0,k+1):
                    summ=summ+(e2-e1)**(k-m) * (e3-e1)**m
                temp=prod*summ
                out=out + temp
                if k>4 and abs(temp) < abs(out)*tol:
                    return out
        if abs(temp) < 10**-100:
            return out
        # appairantly the series didn't converge
        e1,e2,e3 = arrange(e1,e2,e3)
        return ((np.exp(N*e3)/(e3**j)) +                  e2*e1*f2(j+1,e1,e2,N)-e3*f2(j,e1,e2,N))/                 ((e2-e3)*(e3-e1)) 
                
    else: # 1 close to 2 but not to e3          
        return ((np.exp(N*e3)/(e3**j)) +                 e2*e1*f2(j+1,e1,e2,N)-e3*f2(j,e1,e2,N))/                 ((e2-e3)*(e3-e1))

