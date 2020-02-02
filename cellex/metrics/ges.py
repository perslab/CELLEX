import numpy as np
import pandas as pd
import time
import datetime
from .esw_star import esw_star
from ..utils.compute_pvalues import compute_pvalues
from ..summarydata import SummaryData

def _ges(data: pd.DataFrame, mean: pd.DataFrame, nnz: pd.DataFrame, n_cells_per_anno: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    """Computes Gene Enrichment Score ES weights for each gene / cell-type

    Parameters
    ----------
    data : DataFrame
        Expression values per gene / cell.

    mean : DataFrame
        Mean expression per gene / annotation group.

    nnz : DataFrame
        Number of nonzero expression values per gene / annotation group.

    n_cells_per_anno : DataFrame
        Number of cells per annotation group.

    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    enrichment : ndarray
        ES weights
    
    TODO:
        Filter data
        Fix weird computations:
        * np.allclose(np.sum(df.values[:,:],axis=1), (mean_overall*n_cells_total)) >>> True
        * np.allclose(((nnz_overall / n_cells_total) * n_cells_total), nnz_overall) >>> True
    """
    
    # Compute nnz and mean
    mean = mean.values
    nnz = nnz.values

    # number of cells per cluster
    cluster_sizes = n_cells_per_anno.values

    # Non-zeros and mean over all cells
    nnz_overall = np.count_nonzero(data.values,axis=1)
    mean_overall = np.mean(data.values,axis=1)

    # number of cells / columns
    n_cells_total = data.shape[1]

    # Compute fraction of non_zero in clusters and overall
    f_nnz = nnz / cluster_sizes
    f_nnz_overall = nnz_overall / n_cells_total

    # Means and fraction non-zero values in other clusters (per cluster)
    # "((mean_overall * n_cells_total)[None].T - (mean * cluster_sizes))" => means for all genes scaled by n_cells_total minus
    # mean for each cluster scaled by cluster size returns mean of cells outside cluster i.
    # "(n_cells_total - cluster_sizes)" => an array of cell count outside cluster i.
    mean_other = ((mean_overall * n_cells_total)[None].T - (mean * cluster_sizes)) / (n_cells_total - cluster_sizes)
    f_nnz_other = ((f_nnz_overall * n_cells_total)[None].T - (f_nnz * cluster_sizes)) / (n_cells_total - cluster_sizes)

    e1 = 0.01 # org: 0.1
    e2 = 1e-100 # org: 0.01

    enrichment = ((f_nnz + e1) / (f_nnz_other + e1)) * ((mean + e2) / (mean_other + e2))

    # filter values close to 0
    eps_machine = np.finfo(float).eps
    tol = eps_machine**0.5
    enrichment[np.isclose(enrichment, 0, rtol=0, atol=tol)] = 0.

    return enrichment

def ges(stats: SummaryData, verbose: bool=False, compute_meta: bool=False):
    """Compute Gene Enrichment Score

    GES is based on the eponymous metric described in:

        Zeisel, et al. Molecular Architecture of the Mouse Nervous System. 
        Cell 174, 999- 1014.e22 (2018).

    and implemented as MarkerSelection in Cytograph. Code available at:

        github(.)com/linnarsson-lab/cytograph


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
        print("Computing GES ...")

    # df = stats.data
    idx_labels = stats.data.index
    col_labels = stats.mean.columns.values
    key = "ges."

    results = {}

    if verbose:
        print("    esw ...")
    esw = _ges(stats.data, stats.mean, stats.n_nonzero, stats.n_cells_per_anno, verbose)
    esw_df = pd.DataFrame(data=esw, index=idx_labels, columns=col_labels)
    results[(key + "esw")] = esw_df

    if compute_meta:
        esw_null = _ges(stats.data, stats.mean_null, stats.n_nonzero_null, stats.n_cells_per_anno_null, verbose=verbose)
        pvals = compute_pvalues(esw, esw_null, verbose)
        pvals_df = pd.DataFrame(pvals, idx_labels, col_labels)
        results[(key + "esw_null")] = pd.DataFrame(esw_null, idx_labels, col_labels)
        results[(key + "pvals")] = pvals_df
        results[(key + "esw_s")] = esw_star(esw_df, pvals_df, verbose)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))

    return results
