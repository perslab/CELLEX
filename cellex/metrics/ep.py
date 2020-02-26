import numpy as np
import pandas as pd
import time
import datetime
from .esw_star import esw_star
from ..utils.compute_pvalues import compute_pvalues
from ..summarydata import SummaryData

def _ep(mean: pd.DataFrame, verbose: bool=False):
    """Computes Expression Proportion ES weights for each gene / cell-type

    Parameters
    ----------
    mean : DataFrame
        Mean expression per gene / annotation group.

    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    specificity : ndarray
        ES weights

    """

    n_genes = mean.shape[0]

    # n_genes x n_annotations
    expr_mean = mean.values

    # Scale expression values by sum of cell expression:
    # every column should sum to 1.
    # axis=0, i.e. sum column values
    # Shape: n_genes x n_annotations
    expr_mean = expr_mean / np.sum(expr_mean, axis=0)

    # n_genes x 1
    # axis=1, i.e. sum row values. Reshape to ensure vector-shape
    expr_mean_sum = np.sum(expr_mean,axis=1).reshape((n_genes,1))

    # Compute Specificity
    # mean_x / (sum_means_all + eps)
    specificity = expr_mean / (expr_mean_sum + 1e-12)

    return specificity

def ep(stats: SummaryData, verbose: bool=False, compute_meta: bool=False):
    """Compute Expression Proportion

    EP is based on the specificity calculations described in:

        Skene, et al. Genetic identification of brain cell types underlying 
        schizophrenia. Nat. Genet. 50, 825–833 (2018)
    
    and implemented in the EWCE R package. Code available at:
    
        github(.)com/NathanSkene/EWCE


    Parameters
    ----------
    summarydata : SummaryData
        Summary data computed from raw data using specified annotation.

    verbose : bool, optional (default: False)
        Print progress report.
    
    compute_meta : bool, optional (default: False)
        Compute meta results.

    Returns
    -------
    results : dict
        Dictionary containing all computed ESw and meta results, e.g. pvals
    
    """
    
    start = 0

    if verbose:
        start = time.time()
        print("Computing EP ...")

    df = stats.mean
    idx_labels = df.index
    col_labels = stats.mean.columns.values
    key = "ep."

    results = {}

    if verbose:
        print("    esw ...")
    esw = _ep(df, verbose)
    esw_df = pd.DataFrame(esw, idx_labels, col_labels)
    results[(key + "esw")] = esw_df

    if compute_meta:
        esw_null = _ep(stats.mean_null, verbose)
        pvals = compute_pvalues(esw, esw_null, verbose)
        pvals_df = pd.DataFrame(pvals, idx_labels, col_labels)
        results[(key + "esw_null")] = pd.DataFrame(esw_null, idx_labels, col_labels)
        results[(key + "pvals")] = pvals_df
        results[(key + "esw_s")] = esw_star(esw_df, pvals_df, verbose)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))

    return results
