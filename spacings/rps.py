import numpy as np
from collections import namedtuple

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

def rps(x):
    '''
    '''
    ts = rps_ts(x)

    return RPStestResult(ts, 0.0)
