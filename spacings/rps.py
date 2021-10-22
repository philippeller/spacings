import numpy as np
from collections import namedtuple
from scipy import interpolate
from scipy.stats import distributions
from numpy import genfromtxt
import pkg_resources
import warnings

RPStestResult = namedtuple('RPStestResult', ('statistic', 'pvalue'))

degree = int(np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/long_1_RPS_table_fit_smoothing_degree.csv'), delimiter=',', skip_header=1)[1])
data_knots = np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/long_1_RPS_table_fit_Bspline_knots.csv'), delimiter=',', skip_header=1)
data_coefficients = np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/long_1_RPS_table_fit_Bspline_coefficients.csv'), delimiter=',', skip_header=1)
p_vals = np.genfromtxt(pkg_resources.resource_filename('spacings', 'rps_tables/long_1_RPS_table_fit_p_values.csv'), delimiter=',', skip_header=1)

models = []
for i in range(data_knots.shape[0]):
    mask = data_knots[i, 1:] >= 0
    models.append(interpolate.BSpline(data_knots[i, 1:][mask], data_coefficients[i, 1:][mask], degree))

def get_x_vals(N):
    return np.array([f(np.log10(N)) for f in models])

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
    x_vals = get_x_vals(N)
    
    specific_cdf = interpolate.PchipInterpolator(np.concatenate([[0,], x_vals, [1.]]), np.concatenate([[0,], p_vals, [1.]]))

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

    if np.any(np.diff(cdfvals) == 0):
        return RPStestResult(0, 0)

    ts = transformed_rps_ts(cdfvals)

    if N > 1000:
        warnings.warn('p-value calculation only implemented for sample sizes <= 1000')

    p = float(specific_RPS_norm_edge_p_value(N)(ts))

    return RPStestResult(ts, p)