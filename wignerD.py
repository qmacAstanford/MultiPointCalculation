
# coding: utf-8

#   This class allows you to quickly calculate Wigner D matrix values.  I uses the symbolic lookup Wigner D functionality available from sympy.physics.  It will be slow when first calculating a set of m,mp, and l values but will hence-forth become fast and evaluate quickly for different angles.  See the file "Wigner_Convention" for the conventions I use.
# Examples: First you define a set, then use it
# 
# >> import wignerD as wd
# >> myset = wd.wigner_d_vals()
# >> ans = myset.get(l,mp,m,alpha,beta,gamma)
# 
# You can also use a faster algorithm if you have gamma=0.
# 
# >> ans = myset.axial(l,mp,m,alpha,beta)

# In[1]:

from sympy.physics.quantum.spin import Rotation
from sympy import symbols
import numpy as np
from sympy.utilities import lambdify
class wigner_d_vals:
    def __init__(self):
        self.fun = {} #dictionary of functions, key: (l,mp,m)
        self.fun_axial = {} # axial symmetry
        
    def axial(self,l,mp,m,alpha,beta):
        key = (l,mp,m)
        if abs(mp)>l or abs(m)>l:
            raise Exception('m must be <= l')
        if key in self.fun_axial:
            fun=self.fun_axial[key]
        else:
            al, be = symbols('alpha beta')
            fun=lambdify((al,be),Rotation.D(l,mp,m,al,be,0).doit(),"numpy")
            self.fun_axial[key]=fun
        return fun(alpha,beta)    
            
    def get(self,l,mp,m,alpha,beta,gamma):
        key = (l,mp,m)
        if key in self.fun:
            fun=self.fun[key]
        else:
            al, be, ga = symbols('alpha beta gamma')
            fun=lambdify((al,be,ga),Rotation.D(l,mp,m,al,be,ga).doit(),"numpy")
            self.fun[key]=fun
        return fun(alpha,beta,gamma)             
        
    def get_mtrx(self,alpha,beta,gamma,nlam):
        out=[] # list of matricies [l][mp-l,m-l]
        for l in range(0,nlam):
            print('l',l)
            D = np.zeros((2*l+1,2*l+1),dtype='complex')
            for m in range(-l,l+1):
                for mp in range(-l,l+1):
                    # note that I may use m and mp bacward from sympy docs
                    D[mp-l,m-l]=self.get(l,mp,m,alpha,beta,gamma)
            out.append(D)
        return out

