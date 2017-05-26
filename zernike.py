from operator import xor
import numpy as np

def fringe_1d_to_2d(j):
    '''
    Convert a Fringe 1D index to the corresponding
    classical radial, azimuthal (n, m) pair.
    '''
    assert j > 0, 'j must be > 0! (Piston is j=1)'
    
    d = np.floor(np.sqrt(j-1)) + 1
    if not is_even(d**2 - j):
        m = np.ceil( (d**2 - j) / 2. )
    else:
        m = np.ceil( (-d**2 + j -1) / 2. )
    n = 2 * (d - 1) - abs(m)
    return int(n), int(m)

def fringe_2d_to_1d(n,m):
    '''
    Convert an (n, m) pair to the Fringe 1D index.
    
    Note: this assumes (n, m) are the classical
    radial, azimuthal degree pair -- not the
    Wyant (n, m).
    
    '''
    sgn = lambda x : np.sign(x) if x !=0 else 1 # treat 0 as positive
    
    _nm_validation(n,m)
    
    j = ( (n + abs(m)) / 2 + 1)**2 - 2 * abs(m) + ( 1 - sgn(m) ) / 2
    return int(j)

def noll_1d_to_2d(j):
    pass

def noll_2d_to_1d(n,m):
    pass

def classical_nm_to_Wyant(n,m):
    '''
    Convert classical radial, azimuthal degrees
    into Wyant's all-positive degrees. Note that
    this returns an additional parameter to break
    the degeneracy in the Wyant convention.
    '''
    _nm_validation(n,m)
    
    nprime = ( n + abs(m) ) // 2
    sign = np.sign(m)
    if sign > 0:
        azimuthal = 'cos'
    elif sign < 0:
        azimuthal = 'sin'
    else:
        azimuthal = None
    return nprime, abs(m), azimuthal

def noll_to_fringe_normalization(n,m):
    pass

def fringe_to_noll_normalization(n,m):
    pass

def is_even(x):
    return x % 2 != 0

def _nm_validation(n,m):
    assert n >= 0, 'n must be positive!'
    assert abs(m) <= n, 'Absolute value of m cannot exceed n!'
    assert not xor(is_even(n),is_even(m)), '(n,m) must be both even or both odd!'

def _check_fringe(nterms=36):
    '''
    Check that fringe_2d_to_1d and fringe_1d_to_2d
    are consistent.
    '''
    for j in range(1,nterms+1):
        n, m = fringe_1d_to_2d(j)
        foundj = fringe_2d_to_1d(n, m)
        assert j == foundj, 'Fringe indices don\'t match at j={} (found j={})'.format(j,foundj)