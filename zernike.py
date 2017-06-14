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
    '''
    Given a Noll index (j >= 1), return the (n,m) pair.
    
    The treatment here closely follows Robert Gray's ZernikeCalc.m
    '''
    
    assert j > 0, 'j must be > 0! (Piston is j=1)'
    
    # Find n
    n = int(np.ceil( (-1 + np.sqrt(1 + 8 * j)) / 2 - 1 ))
    
    # Find sign of m
    if is_even(j):
        sign = -1
    else:
        sign = 1
        
    # Find m
    k = (n + 1) * (n + 2) / 2.
    m = sign * int(n - 2 * np.floor((k - j) / 2. ) )
    
    return n, m

def noll_2d_to_1d(n,m):
    '''
    Given an (n,m) pair, return the corresponding
    Noll index.
    '''
    _nm_validation(n,m)
    
    j0 = n * (n + 1.) / 2. + 1
    
    if is_even(j0 + n):
        sign = -1
    else:
        sign = 1
        
    # There has to be a more elegant way to do this
    mvals = np.arange(-n,n+1,2)
    morder = np.lexsort((sign * np.sign(mvals), np.abs(mvals)))
    
    return j0 + list(mvals[morder]).index(m)

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

def noll_normalization(n,m):
    if m == 0:
        return np.sqrt(n + 1)
    else:
        return np.sqrt(2 * (n + 1))

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

def _check_noll(nterms=36):
    '''
    Check that noll_2d_to_1d and noll_1d_to_2d
    are consistent.
    '''
    for j in range(1,nterms+1):
        n, m = noll_1d_to_2d(j)
        foundj = noll_2d_to_1d(n, m)
        assert j == foundj, 'Noll indices don\'t match at j={} (found j={})'.format(j,foundj)