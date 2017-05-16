
# coding: utf-8

# Here are some functions that are useful and have some numerical limits that have been fixed.

# In[1]:

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

# In[2]:

def expl(x,n,maxn=15):
    if abs(x)<0.001:
        out=0
        for ii in range(n,maxn):
            out=out+x**ii/np.math.factorial(ii)
    else:
        out=np.exp(x)
        for ii in range(0,n):
            out=out-x**ii/np.math.factorial(ii) 
    return out


# In[3]:

# sequence of integers from first to last
# for example seq(2,5)=2,3,4,5
@jit(nopython=True)
def seq(first,last):
    if first<=last:
        return np.arange(first,last+1)
    else:
        return np.arange(first,last-1,-1)


# In[4]:

#@jit(nopython=True)
def factorial():
    out=int(1)
    fact = np.ones(26)
    for ii in range(1,25+1):
        out=out*ii
        fact[ii] = out
    return fact
fact=factorial()


# 
# \[
# \frac{\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)}{\left(\epsilon_{1}-\epsilon_{2}\right)}
# \]

# In[5]:

@jit(nopython=True)
def f2(j,e1,e2,N):
    tol=0.001
    nmax=10
    if abs(e2-e1)<tol:
        x=(e2-e1)/2
        y=(e2+e1)/2
        if j == 1:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-(-N*y)**m / fact[m]
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-2))
        elif j == 2:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-(n+2-m)*(-N*y)**m / fact[m]
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-3))
        elif j== 3:
            out=0.0
            for n in range(0,nmax+1,2):
                temp=0.0
                for m in seq(0,n+1):
                    temp=temp-0.5*(n+3-m)*(n+2-m)*(-N*y)**m / fact[m]
                out=out+temp*(x**n)*np.exp(N*y)*(y**(-n-4))
    else:
        out=(np.exp(N*e1)/(e1**j) - np.exp(N*e2)/(e2**j))/(e1-e2)
    return out


# In[6]:

@jit(nopython=True)
def arrange(a,b,c):
    ab = abs(a-b)
    bc = abs(b-c)
    ac = abs(a-c)
    if ab <= min(bc,ac):
        return (a,b,c)
    elif ac <= min(ab,bc):
        return (a,c,b)
    else:
        return (b,c,a)
        


# \[
# \frac{\left(\frac{e^{N\epsilon_{3}}}{\epsilon_{3}^{j}}-\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}\right)\epsilon_{1}+\left(\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}-\frac{e^{N\epsilon_{3}}}{\epsilon_{3}^{j}}\right)\epsilon_{2}+\left(\frac{e^{N\epsilon_{2}}}{\epsilon_{2}^{j}}-\frac{e^{N\epsilon_{1}}}{\epsilon_{1}^{j}}\right)\epsilon_{3}}{\left(\epsilon_{1}-\epsilon_{2}\right)\left(\epsilon_{2}-\epsilon_{3}\right)\left(\epsilon_{3}-\epsilon_{1}\right)}
# \]

# In[13]:

@jit(nopython=True)
def f3(j,e1,e2,e3,N):
    e1,e2,e3 = arrange(e1,e2,e3)
    tol=0.01
    if abs(e1-e2) >= tol: # none close to eachother
        out=((np.exp(N*e3)/(e3**j) - np.exp(N*e2)/(e2**j))*e1+             (np.exp(N*e1)/(e1**j) - np.exp(N*e3)/(e3**j))*e2+             (np.exp(N*e2)/(e2**j) - np.exp(N*e1)/(e1**j))*e3)/            ((e1-e2)*(e2-e3)*(e3-e1))
    elif( max(abs(e1-e2),abs(e2-e3),abs(e3-e1)) <= tol):
        # all close to eachother
        if j==1:
            out=-np.exp(N*e1)*((N*e1)**2 -2*N*e1 +2)/(2*e1**3)
        elif j==2:
            out=-np.exp(N*e1)*((N*e1)**2 -4*N*e1 +6)/(2*e1**4)
    else:
        # 1 close to 2 but not to e3
        out=((np.exp(N*e3)/(e3**j)) + e2*e1*f2(j+1,e1,e2,N)-e3*f2(j,e1,e2,N))/            ((e2-e3)*(e3-e1))
    return out

