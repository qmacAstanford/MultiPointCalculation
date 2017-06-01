import numpy as np
import WLCgreen as wlc


# Find all shared Eignevalues
# returns of list of overlaps of the format out[l_2]=l_1
# assumes no repeated values within tol
def IntersectEig(set1,set2,tol=10**-4):
    out=np.zeros(len(set1),dtype='int')-1
    # The following is a O(n^2) precess
    # This can be accelerated useing the fact that they are sorted
    for l1 in range(0,len(set1)):
        for l2 in range(0,len(set2)):
            if abs(set1[l1]-set2[l2])<tol:
                out[l1]=l2
                # check to see if next is closer
                if(l2<len(set2)-1):
                    if abs(set1[l1]-set2[l2+1])<abs(set1[l1]-set2[l2]):
                        out[l1]=l2+1

                break
    

   # Assumes two reverse ordered lists with no repeated values
   # 
   # l1=0
   # l2=0
   # while l1<len(set1) and l2<len(set2):
   #     if abs(set1[l1]-set2[l2])<tol:
   #         out[l1]=l2
   #         l1=l1+1
   #         l2=l2+1
   #         continue  
   #     if set1[l1].real >= set2[l2].real + tol:
   #         l1=l1+1
   #     elif set1[l1].real <= set2[l2].real - tol:
   #         l2=l2+1
   #     elif set1[l1].imag >= set2[l2].imag:
   #         l1=l1+1
   #     else
   #         l2=l2+1
        
    return out

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
        self.eig, self.res, self.a, self.b = self.calc_resData_G()
        self.G0, self.dG0 = self.calcZeroPterms()
        self.otherG = {}
       
    def __str__(self):
        temp=np.get_printoptions()
        np.set_printoptions(precision=3)
        print('---Propagator ',self.name,'---')
        print('nlam',self.nlam)
        print('eig:')
        print(self.eig[0:5])
        print('res of lam=lam0=mu:')
        print(self.res[0:5][sel.mu,self.mu])
        print('a of lam=lam0=mu:')
        print(self.a[0:5][self.mu,self.mu])
        print('b of lam=lam0=mu:')
        print(self.b[0:5][self.mu,self.mu])
        print('Other G:')
        for key, value in self.otherG.items():
            print(key,end='')
        print('')
        np.set_printoptions(**temp)
        return('G(K='+str(self.K)+',mu='+str(self.mu)+')= <'+str(self.name)+'>')
     
    def printOther(self,key=None):
        temp=np.get_printoptions()
        np.set_printoptions(precision=3)
        print('Other Name:',key)
        print('overlapping_other_l[l_self]=l_other')
        print(self.otherG[key][1][0:3])
        print('G(lam=self.mu,lam0=self.mu)=')
        print(self.otherG[key][0][self.mu,self.mu,0:3])
        np.set_printoptions(**temp)
    
    def printOthers(self):
        for key in self.otherG:
            self.printOther(key)
              
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
        a=[]
        b=[]
        
        for l in range(0,ORDEig):
            if l<mu:
                eig.append(np.NaN)
                res.append(np.NaN)
                a.append(np.NaN)
                b.append(np.NaN)
            else:
                resTemp, aTemp, bTemp,\
                        invG, dinvG, ddinvG, dddinvG =\
                        wlc.LurentifyG(K,eigvals[l],l,mu,nlam=10,d=3,\
                                       lamMax=500,\
                                       cutoff=10**-11,lowKcut=10**-7)
                #res.append(wlc.residues(K,eigvals[l],l,mu,nlam=nlam,d=d))
                eig.append(eigvals[l])
                res.append(resTemp)
                a.append(aTemp)
                b.append(bTemp)
        return eig, res, a, b
    
    
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
   
    def get_G_laplace(self,p,lam0,lam,mu):
        mu =self.mu
        K = self.K
        lamMax = self.lamMax
        nlam = self.nlam
        jp=wlc.get_jp(p,mu,K,lamMax)
        jm=wlc.get_jm(p,mu,K,lamMax)
        W=wlc.get_W(p,mu,K,lamMax,jp,jm)
        return wlc.get_G(p,lam,lam0,mu,K,jm,jp,W)

    # Evaluate this G at eigenvalues of another G
    # If they share eigenvalues then retun NaN
    # Output format: [self.lam0,self.lam,other.l]
    def get_G_matrix(self,other,tol=10**-4):
        mu =self.mu
        K = self.K
        lamMax = self.lamMax
        nlam = self.nlam
        neig=len(other.eig)
        
        overlapping_other_l = IntersectEig(self.eig,other.eig,tol)
        
        G_at_other = np.zeros((nlam,nlam,neig),dtype=type(1+1j))*np.NaN  
        
        for ii in range(other.mu,neig):
            if overlapping_other_l[ii]==-1: # if no overlaps
                p=other.eig[ii]
                jp=wlc.get_jp(p,mu,K,lamMax)
                jm=wlc.get_jm(p,mu,K,lamMax)
                W=wlc.get_W(p,mu,K,lamMax,jp,jm)
                for lam0 in range(abs(mu),nlam):
                    for lam in range(abs(mu),nlam):
                        G_at_other[lam0,lam,ii]=wlc.get_G(p,lam,lam0,\
                                                          mu,K,jm,jp,W)
            else: # a.k.a overlapping eigenvalue
                G_at_other[:,:,ii]=np.NaN 
        return G_at_other, overlapping_other_l
    
    def isAnEigenvalue(self,p,tol=10**-4):
        for l in range(self.mu,self.ORDEig):
            if abs(p-self.eig[l])<tol:
                return l
        return -1 # not in list
    
    # Add ability to lookup G at eigenvalues of another G
    def add_another_G(self,name,eig):
        if name in self.otherG:
            print(name,' is already in ',self.name)
        else:
            self.otherG.update({name,self.get_G_matrix(eig)})
            
    # Look up G[lam0,lam] evaluate at other G's lth eigenvalue
    # If you don't already know it, add to dictionary
    # lam0, lam are the angular indices of the self
    # other_l is the eigenvalue index of the other propagator
    def G_other(self,lam0,lam,other,other_l):
        name = other.name        
        if name in self.otherG:
            G_at_other, overlapping_other_l =self.otherG[name]
            return G_at_other[lam0,lam,l], overlapping_other_l
        else:
            self.otherG.update({name: self.get_G_matrix(other)})
            G_at_other, overlapping_other_l =self.otherG[name]
            return G_at_other[lam0,lam,l], overlapping_other_l
    
    # Same as G_other except return residues of all other_l 
    def G_others(self,lam0,lam,other):
        name = other.name
        
        if name in self.otherG:
            G_at_other, overlapping_other_l =self.otherG[name]
            return G_at_other[lam0,lam,:], overlapping_other_l
        else:
            #print('self.K',self.K,'other.K',other.K)
            self.otherG.update({name: self.get_G_matrix(other)})
            G_at_other, overlapping_other_l =self.otherG[name]
            return G_at_other[lam0,lam,:], overlapping_other_l
            


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
            self.mu_dict[abs(mu)]=propagator(name,\
                                             np.linalg.norm(self.kvec),\
                                             abs(mu),\
                                             nlam=self.nlam,\
                                             lamMax=self.lamMax,\
                                             ORDEig=self.ORDEig,\
                                             d=self.d)
            return self.mu_dict[abs(mu)]        



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
            self.prop_dict[key]=vectorPropagator(kapprox,\
                                                 self.nlam,\
                                                 self.lamMax,\
                                                 self.ORDEig,\
                                                 self.d)
            return self.prop_dict[key]
        
        
