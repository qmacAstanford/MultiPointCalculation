import numpy as np
import MultiPoint.WLCgreen as wlc

def IntersectEig(set1,set2,tol=10**-4):
    """
    Find all shared Eignevalues
    assumes no repeated values within tol
    Args:
        set1 (interable): list of values
        set2 (interable): list of values
    Returns:
        out: List of overlaps of the format out[l_2]=l_1
    """
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
    return out

def IntersectEig2(set1,set2,tol=10**-8):
    """
    Find all shared Eignevalues
    Finds values in set1 and set2 that are within tol of eachother.
    Assumes no repeated values within tol
    Args:
        set1 (interable): list of values
        set2 (interable): list of values
    Returns:
        set1only (numpy array): length of set1. 1 if no overlaps for that
        eigenvalue, 0 if there is.

        set2only (numpy array): length of set2. 1 if no overlaps for that
        eigenvalue, 0 if there is.

        paris list((int,int)): returns a list of tuples of indicies of overlap
        (l1,l2)
    """
    set1only=np.ones(len(set1),dtype='int')
    set2only=np.ones(len(set1),dtype='int')
    pairs=[]
    for l1 in range(0,len(set1)):
        for l2 in range(0,len(set2)):
            if set2only[l2] == 0:
                continue
            if abs(set1[l1]-set2[l2]) < tol:
                set1only[l1]=0
                set2only[l2]=0
                pairs.append((l1,l2))
                break 
    return set1only, set2only, pairs
                
def IntersectEig3(set1,set2,set3,tol=10**-10):
    """
    Find all shared Eignevalues
    Finds values in set1, set2, and set3 that are within tol of eachother.
    Assumes no repeated values within tol
    Args:
        set1 (interable): list of values
        set2 (interable): list of values
        set3 (interable): list of values
    Returns:
        set1only (numpy array): length of set1. 1 if no overlaps for that
        eigenvalue, 0 if there is.

        set2only (numpy array): length of set2. 1 if no overlaps for that
        eigenvalue, 0 if there is.

        set3only (numpy array): length of set3. 1 if no overlaps for that
        eigenvalue, 0 if there is.

        paris12 list((int,int)): returns a list of tuples of indicies of
        overlap (l1,l2)

        paris23 list((int,int)): returns a list of tuples of indicies of
        overlap (l2,l3)

        paris31 list((int,int)): returns a list of tuples of indicies of
        overlap (l3,l1)

        triplets list((int,int,int)): list of (l1,l2,l3) that overlap
    """
    set1only=np.ones(len(set1),dtype='int')
    set2only=np.ones(len(set1),dtype='int')
    set3only=np.ones(len(set1),dtype='int')
    pairs12=[]
    pairs23=[]
    pairs31=[]
    triplets=[] 
    for l1 in range(0,len(set1)):
        for l2 in range(0,len(set2)):
            if set2only[l2] == 0:
                continue
            if abs(set1[l1]-set2[l2]) < tol:
                set1only[l1]=0
                set2only[l2]=0
                thirdOverlap=0
                for l3 in range(0,len(set3)):
                    if set3only[l3] == 0:
                        continue
                    if abs(set1[l1]-set3[l3])< tol or abs(set2[l2]-set3[l3])< tol:
                        set3only[l3]=0
                        triplets.append((l1,l2,l3))
                        thirdOverlap=1
                if thirdOverlap==0:
                    pairs12.append((l1,l2))
                break
                
        for l3 in range(0,len(set3)):
            if set3only[l3] == 0:
                continue
            if abs(set1[l1]-set3[l3]) < tol:
                set1only[l1]=0
                set3only[l3]=0
                pairs31.append((l3,l1))
                break
             
    for l2 in range(0,len(set1)):
        if set2only[l2] == 0:
            continue
        for l3 in range(0,len(set3)):
            if set3only[l3] == 0:
                continue
            if abs(set3[l3]-set2[l2]) < tol:
                set3only[l3]=0
                set2only[l2]=0
                pairs23.append((l2,l3))
                break
    return set1only, set2only, set3only, pairs12, pairs23, pairs31, triplets
                              
        
                
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
    def __init__(self,name,K,mu,nlam=10,lamMax=500,ORDEig=25,d=3,
             resCutoff=10**-11,lowKcut=10**-7,eigCutoff=800):
        self.name = name
        self.K = K
        if d != 3 and mu <0:
            print('Error in propagator:')
            print('I assumed that G(mu)=G(-mu).')
            print('This appears to be only true if d=3.')
        self.mu = abs(mu)
        self.lamMax = lamMax
        self.resCutoff=resCutoff
        self.lowKcut=lowKcut
        self.nlam = nlam
        self.ORDEig = ORDEig
        self.d =d
        self.eigCutoff=eigCutoff
        self.eig, self.res, self.a, self.b = self.calc_resData_G()
        self.G0, self.dG0 = self.calcZeroPterms()
        self.otherG = {}
       
    def __str__(self):
        temp=np.get_printoptions()
        np.set_printoptions(precision=3)
        print('---Propagator ',self.name,'---')
        print('nlam',self.nlam)
        print('mu',self.mu)
        print('K',self.K)
        print('eig:')
        print(self.eig[0:5])
        print('res of lam=lam0=mu:')
        for l in range(0,5):
            print(self.res[l][self.mu,self.mu])
        print('a of lam=lam0=mu:')
        for l in range(0,5):
            print(self.a[l][self.mu,self.mu])
        print('b of lam=lam0=mu:')
        for l in range(0,5):
            print(self.b[l][self.mu,self.mu])
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

        eigvals = wlc.wlc_eigvals(K,ORDEig,mu,d=d,cutoff=self.eigCutoff)
        eig=[]
        res=[]
        a=[]
        b=[]
        
        for l in range(0,ORDEig):
            if l<mu:
                eig.append(np.NaN)
                res.append(np.zeros((self.nlam,self.nlam))*np.NaN)
                a.append(np.zeros((self.nlam,self.nlam))*np.NaN)
                b.append(np.zeros((self.nlam,self.nlam))*np.NaN)
            else:
                resTemp, aTemp, bTemp,\
                        invG, dinvG, ddinvG, dddinvG =\
                        wlc.LurentifyG(K,eigvals[l],l,mu,nlam=nlam,d=d,\
                                       lamMax=self.lamMax,\
                                       cutoff=self.resCutoff,\
                                       lowKcut=self.lowKcut)
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

    # Inverse laplace trasform G(p,K) -> G(N,K)
    # Peformed by summing residues
    # G(N) = \sum_l R(eps_l) exp(N*eps_l)
    def get_G(self,N,lam0,lam):
        eig=self.eig
        res = self.res
        out=0.0
        for l in range(self.mu,len(eig)):
            out=out+res[l][lam0,lam]*np.exp(N*eig[l])
        return out
   
    # Return G(K,p)
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
    # Inputs:
    #     other: propagator object at other
    #     otherOnly: array of 0=overlap, 1=no overlap indexed by lother
    #     tol: returns NaN when within tol of an eigenvalue
    # Outputs:
    #     G_at_other[lam0,lam,lother]: 
    #     dG_at_other[lam0,lam,lother]:
    def get_G_at_other(self,other,tol=10**-13):
        mu =self.mu
        K = self.K
        lamMax = self.lamMax
        nlam = self.nlam
        d=self.d
        neig=len(other.eig)
        
        G_at_other = np.zeros((nlam,nlam,neig),dtype=type(1+1j))*np.NaN  
        dG_at_other = np.zeros((nlam,nlam,neig),dtype=type(1+1j))*np.NaN  
        for lother in range(other.mu,neig):
            # check for overlaps [to avoide divide by zero]
            overlap=self.isAnEigenvalue(other.eig[lother],tol)
            if overlap==-1: # if no overlaps 
                p=other.eig[lother]
                jp=wlc.get_jp(p,mu,K,lamMax,d)
                jm=wlc.get_jm(p,mu,K,lamMax,d)
                djp=wlc.get_djp(p, mu, K, lamMax, jp, d)
                djm=wlc.get_djm(p, mu, K, lamMax, jm, d)
                dwp=wlc.get_dw(jp,djp,mu)
                dwm=wlc.get_dw(jm, djm, mu)
                W=wlc.get_W(p,mu,K,lamMax,jp,jm,d)
                dW=wlc.get_dW(p, mu, K, lamMax, dwp, dwm, W, d)
                for lam0 in range(abs(mu),nlam):
                    for lam in range(abs(mu),nlam):
                        G_at_other[lam0,lam,lother]=wlc.get_G(p,lam,lam0,
                                                          mu,K,jm,jp,W,d)
                        dG_at_other[lam0,lam,lother]=wlc.get_dG(p, lam, lam0,
                                                                mu, K, jm, jp,
                                                                W, dW, dwm,
                                                                dwp,d)
            else: # a.k.a overlapping eigenvalue
                G_at_other[:,:,lother]=np.NaN 
                dG_at_other[:,:,lother]=np.NaN 
        return G_at_other, dG_at_other
    
    def isAnEigenvalue(self,p,tol=10**-4):
        for l in range(self.mu,self.ORDEig):
            if abs(p-self.eig[l])<tol:
                return l
        return -1 # not in list
            
    # Look up G[lam0,lam] evaluate at other G's other_l eigenvalue
    # If you don't already know it, add to dictionary
    # lam0, lam are the angular indices of the self
    # other_l is the eigenvalue index of the other propagator
    # If other_l=None return residues of all other_l
    def G_others(self,lam0,lam,other,other_l=None):
        name = other.name        
        if name not in self.otherG:
            self.otherG.update({name: self.get_G_at_other(other)})
        G_at_other = self.otherG[name][0] # zero for G_at_other

        if other_l==None:
            return G_at_other[lam0,lam,:]  # Return at all other_l
        else:
            return G_at_other[lam0,lam,otherl]

    def dG_others(self,lam0,lam,other,other_l=None):
        name = other.name        
        if name not in self.otherG:
            self.otherG.update({name: self.get_G_at_other(other)})
        dG_at_other = self.otherG[name][1] # 1 for dG_at_other

        if other_l==None:
            return dG_at_other[lam0,lam,:]  # Return at all other_l
        else:
            return dG_at_other[lam0,lam,otherl]
            

class vectorPropagator:
    
    def __init__(self,kvec,nlam=10,lamMax=500,ORDEig=25,d=3):
        self.kvec = kvec
        self.lamMax = lamMax
        self.nlam = nlam
        self.ORDEig = ORDEig
        self.d =d  
        
        self.mu_dict = {} # dictionary of scalar propagators. key: abs(mu)
        self.D= None # Wigner D matrix
        
    def prop(self,mu,precisionPlacesSameProp=8):
        if abs(mu) in self.mu_dict:
            # Already exists, just return it
            return self.mu_dict[abs(mu)]
        else:
            if self.d != 3 and mu <0:
                print('Error in vectorPropagator')
                print('I assumed that G(mu)=G(-mu).')
                print('This appears to be only true if d=3.')

            # name propagator by mu and |k|
            temp=np.get_printoptions()
            np.set_printoptions(precision=precisionPlacesSameProp)
            name=(abs(mu),str(np.linalg.norm(self.kvec)))
            np.set_printoptions(**temp)

            # Calculate
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
        
        
