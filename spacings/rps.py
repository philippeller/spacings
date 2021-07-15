import numpy as np
from collections import namedtuple
from scipy import interpolate
from scipy.stats import distributions
from numpy import genfromtxt
import pkg_resources

RPStestResult = namedtuple('RPStestResult', ('statistic', 'pvalue'))

fit_data_x = np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/RPS_norm_edge_table_fit_spline__p_to_x.csv'), delimiter=',')
fit_data_p = np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/RPS_norm_edge_table_fit_spline__x_to_p.csv'), delimiter=',')

N_fit = fit_data_x[0,1:]

models_N_x = [interpolate.PchipInterpolator(np.log10(N_fit), fit_data_x[i, 1:]) for i in range(1, fit_data_x.shape[0])]
models_N_p = [interpolate.PchipInterpolator(np.log10(N_fit), fit_data_p[i, 1:]) for i in range(1, fit_data_p.shape[0])]

def get_x_vals(N):
    return np.array([f(np.log10(N)) for f in models_N_x])
def get_p_vals(N):
    return np.array([f(np.log10(N)) for f in models_N_p])

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

def transformed_rps_ts(x):
    '''
    Calculate the transformed RPS test statistic

    x : array, sorted
    '''
    ts = 0.
    min_ts = 0.
    while len(x) > 2:
        x = (x - x[0])/(x[-1] - x[0])
        d = np.diff(x)
        ts -= np.sum(np.log(d))
        min_ts -= len(d)*np.log(len(d))
        x = (x[1:] + x[:-1])
        
    return -min_ts/ts

def specific_RPS_norm_edge_p_value(N):
    x2 = fit_data_p[1:,0]
    p1 = fit_data_x[1:,0]
    x1 = get_x_vals(N)
    p2 = get_p_vals(N)
    
    x = np.concatenate([[0], x1, x2, [1]])
    p = np.concatenate([[0], p1, p2, [1]])
    
    idxs = np.argsort(x)
    
    specific_cdf = interpolate.PchipInterpolator(x[idxs], p[idxs])

    return specific_cdf

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

    ts = transformed_rps_ts(cdfvals)

    p = float(specific_RPS_norm_edge_p_value(N)(ts))

    return RPStestResult(ts, p)
