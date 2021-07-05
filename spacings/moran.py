import numpy as np
from collections import namedtuple
from scipy import stats
from scipy.stats import distributions

MorantestResult = namedtuple('MorantestResult', ('statistic', 'pvalue'))

def moran_params(n):
    mu = n * (np.log(n) + np.euler_gamma) - 0.5 - 1/(12*n)
    var = n * (np.pi**2/6 - 1) - 0.5 - 1/(6*n)

    c1 = mu - np.sqrt(var*n/2)
    c2 = np.sqrt(var/(2*n))
    
    return c1, c2, mu, var

def moran(x, cdf='uniform', args=()):
    '''
    Calculates the Moran test and p-value approximation, as described in:
    .. R. C. H. Cheng and M. A. Stephens, “A goodness-of-fit 
       test using moran’s statistic with estimated parameters,
       ”Biometrika, vol. 76, no. 2, pp. 385–392, 1989.

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
        Moran test statistic
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

    ts = -np.sum(np.log(np.diff(cdfvals)))
    c1, c2, mu, var = moran_params(N+1)
    a = (ts - c1)/c2
    p = 1 - stats.chi2(df=N+1).cdf(a)

    return MorantestResult(ts, p)
