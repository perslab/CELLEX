import numpy as np
import pandas as pd

def esw_star(esw: pd.DataFrame, pvals: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    """Compute ESw_star
    
    Computes rank-normalized ESw, i.e. ESw*

    Parameters
    ----------
    esw : DataFrame
        ESw's.
    
    pvals : DataFrame
        Empirical p-values.
    
    Returns
    -------
    esw_ranknorm : DataFrame
        Rank normalized ESw's between 0 and 1.
    """
    if verbose:
        print("    esw_s ...")
    
    # mask matrix of size n_cells x n_genes
    pval_mask = (pvals <= 0.05).values
    binzero_mask = (esw > 0).values
    # invert mask for values to overwrite with NaN
    mask = ~(pval_mask & binzero_mask)
    
    # masked esw matrix of shape n_cells x n_genes
    # significant ESw's are kept, while others are set to NaN
    df_nominal = esw.mask(mask) 

    # rank remaining significant ESw's per column and set NaN's to 0
    df_ranked = (df_nominal.rank(axis=0, method="average")).fillna(value=0)
    
    # divide by number of nominally significant genes group / cell type 
    # to get values between 0 and 1.
    esw_ranknorm = (df_ranked / df_ranked.max())

    return esw_ranknorm
