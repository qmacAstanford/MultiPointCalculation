import numpy as np
from MultiPoint.WLCgreen import wlc_eigvals, residues

def arbitrary_int(L_poly, k_vals, diag, ORDEig = 20, verbose = False):
    """ 
    """
    if hasattr(k_vals,'__len__'):
        return [arbitrary_int(L_poly, kk, diag, ORDEig=ORDEig) for kk in k_vals]
    if verbose:
        print(k_vals)
    mu = 0
    eps = wlc_eigvals(k_vals, ORDEig, mu)
    Res = np.zeros(ORDEig)
    for ell in range(ORDEig):
        Res[ell] = np.real(residues(k_vals, eps[ell], ell, mu, nlam=1)[0,0])

    if verbose:
        print("done with eig/res")
    
    dels = L_poly/float(len(diag)) # path lengh of one bead
    Delta_s_vals = np.arange(len(diag))*dels

    total = 0.0
    for ell in range(len(Res)):
        total += Res[ell] * np.dot(np.exp(eps[ell]*Delta_s_vals), diag)
    if verbose:
        print('done')
    return total*(dels**2)

def test_diagonal_int():
    Ns=300000
    sequence1 = np.random.random(Ns)
    sequence2 = np.random.random(Ns)
    return diagonal_int(sequence1, sequence2)


def diagonal_int(sequence1,sequence2):
    """

    Args:
        sequence1 (array): First sequence
    """
    Ns=len(sequence1)
    assert len(sequence2)==Ns, "sequence1 and sequence2 must have the same length."
    out = np.zeros(Ns)
    for Delta_s in range(0,Ns):
        out[Delta_s] = np.dot(sequence1[0:Ns-Delta_s],sequence2[Delta_s:Ns])
    return out


my_diags = {}

def arb_int(sequence_file, L_poly, k_vals, select='AA'):
    """Numerically intigrates the following (Euler's method)
    .. math::
        \sum_{l}Res_{l}\int_{0}^{N}ds_{2}\int_{0}^{s_{2}}ds_{1}e^{\\epsilon_{l}\left(s_{2}-s_{1}\\right)}\delta_{\alpha_1,\sigma(s_1)}\delta_{\alpha_2,\sigma(s_2)}

    Args:
        sequence_file (string): Filename of sequence of ASCII 1's and 0's
        L_poly (float): Length of polymer in kuhn lengths
        k_vals (float): magnitude of wave vector in inverse kuhn lengths
        select (string): One of 'AA', 'AB', 'BA', 'BB'
    """
    if (sequence_file,select) in my_diags:
        diag = my_diags[(sequence_file, select)]
    else:
        sequence = np.loadtxt(sequence_file)
        if select == 'AA':
            diag = diagonal_int(sequence,sequence)
        elif select == 'AB':
            diag = diagonal_int(sequence,(1.0-sequence))
        elif select == 'BA':
            diag = diagonal_int((1.0-sequence),sequence)
        elif select == 'BB':
            diag = diagonal_int((1.0-sequence),(1.0-sequence))
        my_diags[(sequence_file, select)] = diag

    return arbitrary_int(L_poly, k_vals, diag)
    

def test_arb_int(filename = '../meth'):
    L_poly = 1.7*(10**5)
    
    k_vals = np.logspace(-3,0.2,20)

    out = arb_int(filename, L_poly, k_vals, select = 'AA')

    import matplotlib.pyplot as plt
    plt.semilogx(k_vals,out)

