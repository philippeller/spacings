import numpy as np
from collections import namedtuple
from scipy.stats import distributions

RPStestResult = namedtuple('RPStestResult', ('statistic', 'pvalue'))

def rps_ts(x):
    '''
    Calculate the RPS test statistic

    x : array, sorted
    '''
    ts = 0.    
    while len(x) > 2:
        x = (x - x[0])/(x[-1] - x[0])
        d = np.diff(x)
        ts -= np.sum(np.log(d*len(d)))        
        x = (x[1:] + x[:-1])
    return ts

def rps(x, cdf='uniform', args=()):
    '''
    Calculates the recursive product of spacings (RPS) test

    Parameters:
    -----------
    x : array
        1-D array of observations of random variables
    cdf : string or callable
        If a callable, that callable is used to calculate the cdf.
        If a string, it should be the name of a distribution in `scipy.stats`,
        which will be used as the cdf function.
    args : tuple, sequence, optional
        Distribution parameters for the `cdf`

    Returns:
    --------
    statistic : float
        RPS test statistic
    pvalue :  float
        p-value.

    Examples:
    ---------
    ToDo

    '''

    if isinstance(cdf, str):
        cdf = getattr(distributions, cdf).cdf

    if np.ma.is_masked(x):
        x = x.compressed()

    N = len(x)

    cdfvals = cdf(x, *args)
    cdfvals = np.sort(cdfvals)
    cdfvals = np.concatenate([[0], cdfvals, [1]])

    ts = rps_ts(cdfvals)

    p = 0.0

    return RPStestResult(ts, p)
