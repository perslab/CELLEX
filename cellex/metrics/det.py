import numpy as np
import pandas as pd
import time
import datetime
from .esw_star import esw_star
from ..utils.compute_pvalues import compute_pvalues
from ..summarydata import SummaryData

def _det(mean: pd.DataFrame, var: pd.DataFrame, n_cells: pd.DataFrame, verbose: bool=False):
    """Computes Differential Expression T-statistic ES weights for each gene / cell-type

    Parameters
    ----------
    mean : DataFrame
        Mean expression per gene / annotation group.
    
    var : DataFrame
        Expression variance per gene / annotation group.

    n_cells : DataFrame
        Number of cells per annotation group.

    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    result : ndarray
        ES weights
    
    TODO:
        * Consider replacing args with SummaryStats. Probably not feasible,
            since only tstat uses var and ncells.
        * Consider the consequences of using fillna
        * I think some details are missing in the computations, 
            e.g. subtracting 1 from denominator in var_pooled
        * The way s_p is computed seems off
    """

    ### Compute pooled variance
    variance = var
    variance.fillna(0, inplace=True)
    #n_cells = df.groupby(grouping, axis=1).count() # assuming this is an n_genes x n_annotations
    n_cells = np.array([n_cells.values] * mean.shape[0]) # faster than count
    # n_cells_sum = np.sum(n_cells, axis=1)
    # s_p^2 = sum((sample_size_i - 1) * variance_i) / sum(sample_size_i - 1)
    sd_pooled = np.sqrt(np.sum(((n_cells - 1) * variance.values[:,:]), axis=1) / (np.sum(n_cells - 1, axis=1)))
    
    ### Compute tstat per column
    # should be possible to vectorize using helper function
    # some of these things could be precomputed for the whole matrix
    result = np.zeros(shape=(mean.shape[0], mean.shape[1])) # n_genes, n_cells
    
    for col in mean:
        
        ### Compute X_1 and X_2
        # Separate X and X_other. Use mean values.
        # X is already mean
        # X_other is processed: (mean_other * n_cells_other) / n_cells_sum
        i = mean.columns.get_loc(col)
        
        X_1 = mean.values[:,i]
        
        mean_others = np.delete(mean.values[:,:], i, axis=1)
        
        n_cells_other = np.delete(n_cells, i, axis=1)
        
        X_others = np.sum((mean_others * n_cells_other), axis=1) / (np.sum(n_cells, axis=1))      
        
        ### Compute n_1 and n_other
        # n_1 is n_cells[1]
        # n_other is sum(n_cells[other])
        n_1 = n_cells[:,i]
        
        n_other = np.sum(n_cells_other, axis=1)
        
        ### Compute scaling factor
        sd_pooled_factor = np.sqrt((1/n_1) + (1/n_other))
        
        ### Compute tstat and save result
        tstat = (X_1 - X_others) / (sd_pooled_factor * sd_pooled)
        
        result[:,i] = tstat
    
    return result

def det(stats: SummaryData, verbose: bool=False, compute_meta: bool=False):
    """Compute Differential Expression T-statistic

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
        print("Computing DET ... ")

    idx_labels = stats.data.index
    col_labels = stats.mean.columns.values
    key = "det."

    results = {}

    if verbose:
        print("    esw ...")
    esw = _det(stats.mean, stats.variance, stats.n_cells_per_anno, verbose)
    esw_df = pd.DataFrame(esw, idx_labels, col_labels)
    results[(key + "esw")] = esw_df

    if compute_meta:
        esw_null = _det(stats.mean_null, stats.variance_null, stats.n_cells_per_anno_null, verbose)
        pvals = compute_pvalues(esw, esw_null, verbose)
        pvals_df = pd.DataFrame(pvals, idx_labels, col_labels)
        results[(key + "esw_null")] = pd.DataFrame(esw_null, idx_labels, col_labels)
        results[(key + "pvals")] = pvals_df
        results[(key + "esw_s")] = esw_star(esw_df, pvals_df, verbose)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))

    return results
