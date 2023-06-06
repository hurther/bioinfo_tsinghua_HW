import numpy as np
from scipy.special import logsumexp

def log_mask_zero(a):
    """
    Compute the log of input probabilities masking divide by zero in log.
    """
    a = np.asarray(a)
    with np.errstate(divide="ignore"):
        return np.log(a)
    
    
def log_normalize(a, axis=None):
    """
    Normalize the log probability.
    
    `a`: log probability
    
    `axis`: axis to normalize
    
    return: normalized log probability
    """
    with np.errstate(under='ignore'):
        a_lse = logsumexp(a, axis=axis, keepdims=True)
    a -= a_lse
    
def normalize(a, axis=None):
    """
    Normalize the probability.
    
    `a`: probability
    
    `axis`: axis to normalize
    
    return: normalized probability
    """
    a_sum = a.sum(axis=axis)
    if axis and a.ndim > 1:
        # Make sure we don't divide by zero.
        a_sum[a_sum==0] = 1
        shape = list(a.shape)
        shape[axis] = 1
        a_sum.shape = shape
        
    a /= a_sum
    
def logaddexp(a, b):
    """
    Compute log(exp(a) + exp(b)) in a numerically stable way.
    
    `a`: log probability
    
    `b`: log probability
    
    return: log(exp(a) + exp(b))
    """
    if np.isinf(a) and a < 0:
        return b
    elif np.isinf(b) and b < 0:
        return a
    else:
        return max(a, b) + np.log1p(np.exp(-np.abs(a-b))) 