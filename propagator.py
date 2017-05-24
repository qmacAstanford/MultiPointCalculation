
# coding: utf-8

# In[1]:

import numpy as np
import WLCgreen as wlc


# In[2]:

class propagator:
    
    # Calculate eigenvalues and residues on creation
    def __init__(self,name,K,mu,nlam=10,lamMax=500,ORDEig=25,d=3):
        self.name = name
        self.K = K
        if d != 3 and mu <0:
            print('Error in propagator:')
            print('I assumed that G(mu)=G(-mu).')
            print('This appears to be only true if d=3.')
        self.mu = abs(mu)
        self.lamMax = lamMax
        self.nlam = nlam
        self.ORDEig = ORDEig
        self.d =d
        self.eig, self.res = self.calc_resData_G()
        self.G0, self.dG0 = self.calcZeroPterms()
        self.otherG = {}
        
    # calculate and store eigenvalues and residues
    def calc_resData_G(self):
        mu=self.mu
        ORDEig=self.ORDEig
        nlam=self.nlam
        d=self.d
        K=self.K

        eigvals = wlc.wlc_eigvals(K,ORDEig,mu,d=d)
        eig=[]
        res=[]
        
        for l in range(0,ORDEig):
            if l<mu:
                eig.append(np.NaN)
                res.append(np.NaN)
            else:
                eig.append(eigvals[l])
                res.append(wlc.residues(K,eigvals[l],l,mu,nlam=nlam,d=d))
        return eig, res
    
    # Calculate G(0) and dG(0)/dp
    def calcZeroPterms(self):
        d=self.d
        lamMax=self.lamMax
        nlam=self.nlam
        K=self.K
        mu=self.mu
 
        if K<10**-10: # special K=0 case
            G0=np.zeros((nlam,nlam),dtype=type(1+1j))
            dG0=np.zeros((nlam,nlam),dtype=type(1+1j))
            for lam0 in range(abs(mu),nlam):
                lam=lam0 # zero unless lam=lam0
                # G= 1/( p + lam*(lam+d-2) )
                if lam0==0: # pole at G(0)
                    G0[lam0,lam]=np.inf
                    dG0[lam0,lam]=np.inf
                else:
                    G0[lam0,lam]=1.0/(lam*(lam+d-2))
                    dG0[lam0,lam]=1.0/((lam*(lam+d-2))**2)  
            return G0, dG0
        
        jp=wlc.get_jp(0.0,mu,K,lamMax,d)
        jm= wlc.get_jm_zero(mu,K,lamMax,d)

        djp=wlc.get_djp(0.0,mu,K,lamMax,jp,d)
        djm=wlc.get_djm_zero(mu,K,lamMax,jm,d)


        dwm=wlc.get_dwm_zero(jm,djm,K,lamMax,mu,d)
        dwp=-djp/(jp**2) # This line prevents jit ???

        W=wlc.get_W_zero(mu,K,lamMax,jp,jm,d)
        dW=wlc.get_dW_zero(mu,K,lamMax,dwp,dwm,W,d)
        
        G0=np.zeros((nlam,nlam),dtype=type(1+1j))*np.NaN
        dG0=np.zeros((nlam,nlam),dtype=type(1+1j))*np.NaN
        for lam0 in range(abs(mu),nlam):
            for lam in range(abs(mu),nlam):
                a=1
                G0[lam0,lam]=wlc.get_G_zero(lam,lam0,mu,K,jm,jp,W,d)
                dG0[lam0,lam]=wlc.get_dG_zero(lam,lam0,mu,K,jm,jp,W,dW,dwm,dwp,d)
        return G0, dG0
    
    # Evaluate this G at eigenvalue of different G
    # eig is a list of eigenvalues from another G
    def get_G_matrix(self,eig):
        mu =self.mu
        K = self.K
        lamMax = self.lamMax
        nlam = self.nlam
        neig=len(eig)
        out = np.zeros((nlam,nlam,neig),dtype=type(1+1j))*np.NaN
        for ii in range(0,neig):
            p=eig[ii]
            jp=wlc.get_jp(p,mu,K,lamMax)
            jm=wlc.get_jm(p,mu,K,lamMax)
            W=wlc.get_W(p,mu,K,lamMax,jp,jm)
            for lam0 in range(abs(mu),nlam):
                for lam in range(abs(mu),nlam):
                    out[lam0,lam,ii]=wlc.get_G(p,lam,lam0,mu,K,jm,jp,W)
        return out
    
    #Inverse laplace trasform G(p,K) -> G(N,K)
    # Peformed by summing residues
    # G(N) = \sum_l R(eps_l) exp(N*eps_l)
    def get_G(self,N,lam0,lam):
        eig=self.eig
        res = self.res
        out=0.0
        for l in range(0,len(eig)):
            out=out+res[l][lam0,lam]*np.exp(N*eig[l])
        return out
    
    # Add ability to lookup G at eigenvalues of another G
    def add_another_G(self,name,eig):
        if name in self.otherG:
            print(name,' is already in ',self.name)
        else:
            self.otherG.update({name,self.get_G_matrix(eig)})
            
    # Look up G[lam0,lam] evaluate at other G's lth eigenvalue
    # If you don't already know it, add to dictionary
    def G_other(self,lam0,lam,other,l):
        name = other.name
        eig = other.eig
        
        if name in self.otherG:
            return self.otherG[name][lam0,lam,l]
        else:
            self.otherG.update({name: self.get_G_matrix(eig)})
            return self.otherG[name][lam0,lam,l]
    
    # Same as G_other except return residues of all l 
    def G_others(self,lam0,lam,other):
        name = other.name
        eig = other.eig
        
        if name in self.otherG:
            return self.otherG[name][lam0,lam,:]
        else:
            self.otherG.update({name: self.get_G_matrix(eig)})
            return self.otherG[name][lam0,lam,:]
            


# In[3]:

class vectorPropagator:
    
    def __init__(self,kvec,nlam=10,lamMax=500,ORDEig=25,d=3):
        self.kvec = kvec
        self.lamMax = lamMax
        self.nlam = nlam
        self.ORDEig = ORDEig
        self.d =d  
        
        self.mu_dict = {} # dictionary of scalar propagators
        self.D= None # Wigner D matrix
        
    def prop(self,mu):
        if abs(mu) in self.mu_dict:
            # Already exists, just return it
            return self.mu_dict[abs(mu)]
        else:
            if self.d != 3 and mu <0:
                print('Error in vectorPropagator')
                print('I assumed that G(mu)=G(-mu).')
                print('This appears to be only true if d=3.')
            # Calculate
            name=(abs(mu),self.kvec)
            self.mu_dict[abs(mu)]=propagator(name,                                        np.linalg.norm(self.kvec),                                        abs(mu),                                        nlam=self.nlam,                                        lamMax=self.lamMax,                                        ORDEig=self.ORDEig,                                        d=self.d)
            return self.mu_dict[abs(mu)]        


# In[5]:

# Set of propagators being used in a particular problem.
# This reuses propagators rather than recalculate them.
class prop_set:
    
    def __init__(self,ktol=10**-14,nlam=10,lamMax=500,ORDEig=25,d=3):
        self.ktol = ktol
        self.ndecimals = int(np.ceil(-np.log10(ktol)))
        self.lamMax = lamMax
        self.nlam = nlam
        self.ORDEig = ORDEig
        self.d =d         
        self.prop_dict = {}
    
    def get_vec_prop(self,kvec):
        kapprox  = np.round(np.array(kvec),self.ndecimals)
        key=(kapprox[0],kapprox[1],kapprox[2])
        if key in self.prop_dict:
            return self.prop_dict[key]
        else:
            self.prop_dict[key]=vectorPropagator(kapprox,                                                     self.nlam,                                                     self.lamMax,                                                     self.ORDEig,                                                     self.d)
            return self.prop_dict[key]
        
        


# In[ ]:



