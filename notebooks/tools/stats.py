import numpy as np
from scipy.stats import norm


def fdr( zscores, q=.1, cV=1, invert_zscores=False, mask=None ):
    """
    Adapted from https://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/
    using a default value of cV
    """
    
    if mask is None:
        mask = np.ones(zscores.shape, dtype=bool)

    inv = -1 if invert_zscores else 1
    zscores = inv * zscores
    
    mask *= (zscores != 0)
    zscores = zscores[mask]
    pvals = norm.sf(zscores)
    
    oidx = np.argsort( pvals )
    pvals = pvals[oidx]
    
    V = pvals.size
    idx = np.arange(1, V+1)
    thrline = idx * q / ( V * cV )
    
    select = pvals <= thrline
    if len(pvals[select]):
        thr = np.max(pvals[select])
        zthr = zscores[oidx][select][-1] * inv
    else:
        thr = None
        zthr = None
        
    pcor = pvals * V * cV / idx
    oidx_r = np.argsort(oidx)
        
    padj = np.zeros(len(pvals))
    prev = 1
    
    for i in idx[::-1]:
        padj[i-1] = np.min( [prev, pvals[i-1] * V * cV / i] )
        prev = padj[i-1]
        
    return thr, zthr, pvals, thrline, pcor, padj