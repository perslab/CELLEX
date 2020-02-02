import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy as sp
import time
import datetime
from scipy import stats
from .esw_star import esw_star
from ..utils.compute_pvalues import compute_pvalues
from ..summarydata import SummaryData

def _nsi(mean: pd.DataFrame, verbose: bool=False):
    """Computes Normalized Specificity Index ES weights for each gene / cell-type

    Parameters
    ----------
    mean : DataFrame
        Mean expression per gene / annotation group.

    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    result : ndarray
        ES weights
    
    TODO:
        Filter data
    """
    eps = 0.000000000001

    # number of genes, cells / columns, unique clusters
    n_genes = mean.shape[0]
    n_cells = mean.shape[1]

    # Compute SI
    # Create output matrice for results
    result = np.zeros(shape=(n_genes, n_cells))

    view = mean.values + eps
    selection = np.arange(n_cells)
        
    ### Compute nsi per column
    for j in range(n_cells): # cycle cells
        # Get indexes of other cells
        other_cells = np.delete(selection, j)      
        
        # compute fold change for cell j.
        # N.B. epsilon is added before loop.
        fc = (view[:,[j]]/view[:,other_cells])

        # filter fc values close to 0 and 1, since these are 
        # assumed to be artefacts of adding epsilon
        eps_machine = np.finfo(float).eps
        tol = eps_machine**0.5
        fc[np.isclose(fc, 0, rtol=0, atol=tol)] = 0.
        fc[np.isclose(fc, 1, rtol=0, atol=tol)] = 0.
        
        # compute gene rank for each fold change result, i.e. rank of gene expr per cell
        # we transspose the matrix once, so that we can iterate over cols as rows.
        # we transpose the matrix once more to get back the original cols.
        fc_ranked = np.array([sp.stats.rankdata(col, method="min") for col in fc.T]).T
        
        # normalize: divide each column by number of genes. Subtract 1 in denominator and numerator
        fc_ranked_norm = (fc_ranked - 1) / (fc_ranked.shape[0] - 1) # i.e. rank / n_rows

        # compute average rank based on fold change ranks, i.e. compute mean per row / gene
        fc_mean = fc_ranked_norm.mean(axis=1)
        
        # store gene_rank results (a column) for cell j
        result[:,j] = fc_mean

    return result

def nsi(stats: SummaryData, verbose: bool=False, compute_meta: bool=False):
    """Compute Normalized Specificity Index

    NSI is based on Specificity Index described in:

        Dougherty, et al. Analytical approaches to RNA profiling data for the 
        identification of genes enriched in specific cells. Nucleic Acids Res. 
        38, 4218–4230 (2010).

    and implemented in R. Code available at:

        bactrap(.)org


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
        print("Computing NSI ...")

    df = stats.mean
    idx_labels = df.index
    col_labels = stats.mean.columns.values
    key = "nsi."

    results = {}

    if verbose:
        print("    esw ...")
    esw = _nsi(df, verbose)
    esw_df = pd.DataFrame(esw, idx_labels, col_labels)
    results[(key + "esw")] = esw_df

    if compute_meta:
        esw_null = _nsi(stats.mean_null, verbose)
        pvals = compute_pvalues(esw, esw_null, verbose)
        pvals_df = pd.DataFrame(pvals, idx_labels, col_labels)
        results[(key + "esw_null")] = pd.DataFrame(esw_null, idx_labels, col_labels)
        results[(key + "pvals")] = pvals_df
        results[(key + "esw_s")] = esw_star(esw_df, pvals_df, verbose)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))

    return results
