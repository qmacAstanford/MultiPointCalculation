from sympy.physics.quantum.spin import Rotation
from sympy import symbols
import numpy as np
from sympy.utilities import lambdify
class wigner_d_vals:
    """
    This class allows you to quickly calculate Wigner D matrix values.  I uses the symbolic lookup Wigner D functionality available from sympy.physics.  It will be slow when first calculating a set of m,mp, and l values but will hence-forth become fast and evaluate quickly for different angles.  See the file "Wigner_Convention" for the conventions I use.


    Examples: First you define a set, then use it
    
    .. code-block:: python
        >>> import MultiPoint.wignerD as wd
        >>> myset = wd.wigner_d_vals()
        >>> l=2; m=1; mp=-1; alpha=0.123; beta=0.321; gamma=0.312
        >>> ans = myset.get(l,mp,m,alpha,beta,gamma)
        >>> print(ans)
        (0.0726923757602399-0.013904819324251715j)
        >>> # You can also use a faster algorithm if you have gamma=0.0
        >>> ans = myset.axial(l,mp,m,alpha,beta)
        >>> print(ans)
        (0.07345116118627482+0.00908033118866673j)
    """
    def __init__(self):
        self.fun = {} #dictionary of functions, key: (l,mp,m)
        self.fun_axial = {} # axial symmetry
        
    def axial(self,l,mp,m,alpha,beta):
        """WignerD element

        Args:
            l (int): Angular momentum quatum number
            mp (int): :math:`m'`, z-compoennet in primed system
            mp (int): :math:`m`, z-compoennet in primed system
            alpha (float): :math:`\\alpha` First (z axis) Euler angle
            beta (float): :math:`\\beta` Second (y axis) Euler angle
        """
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
        """WignerD element 

        Args:
            l (int): Angular momentum quatum number
            mp (int): :math:`m'`, z-compoennet in primed system
            mp (int): :math:`m`, z-compoennet in primed system
            alpha (float): :math:`\\alpha` First (z axis) Euler angle
            beta (float): :math:`\\beta` Second (y axis) Euler angle
            gammaa (float): :math:`\\gamma` Third (z axis) Euler angle
        """
        key = (l,mp,m)
        if key in self.fun:
            fun=self.fun[key]
        else:
            al, be, ga = symbols('alpha beta gamma')
            fun=lambdify((al,be,ga),Rotation.D(l,mp,m,al,be,ga).doit(),"numpy")
            self.fun[key]=fun
        return fun(alpha,beta,gamma)             
        
    def get_mtrx(self,alpha,beta,gamma,nlam):
        """WignerD Matrix

        Args:
            alpha (float): :math:`\\alpha` First (z axis) Euler angle
            beta (float): :math:`\\beta` Second (y axis) Euler angle
            gammaa (float): :math:`\\gamma` Third (z axis) Euler angle
            nlam (int): Number of angular quatum numbers, :math:`l`, to include
        """
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

