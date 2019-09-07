import numpy as np
import pandas as pd

def compute_pvalues(scores: np.ndarray, null_scores: np.ndarray, verbose: bool=False) -> pd.DataFrame:
    """
    Compute empirical p-values.
    
    Args:
        scores (ndarray)      : The ES scores computed using actual annotation.
        null_scores (ndarray) : The ES scores computed using permuted annotation.
    
    Returns:
        ndarray of empirical p-values
    """
    if verbose:
        print("    empirical p-values ...")

    #idx_labels = scores.index
    #col_labels = scores.columns.values

    # Sort values of each column
    # copy to prevent side effects
    null_scores_sorted = null_scores.copy()
    null_scores_sorted.sort(axis=0)

    # 1) Takes a pair of matching columns from scores and null_scores_sorted
    # 2) Applies searchsorted on these columns, i.e. for each value in
    #    observation col: at what index should it be placed in the 
    #    null_scores_sorted col to be correctly sorted.
    #    N.B. this index value is equal to the number of null_score vals
    #    less than the observation.
    # N.B. side="right" to match R findInterval behavior
    scores_smaller = np.array([np.searchsorted(null,obs, side="right") for \
                               null,obs in \
                               np.stack((null_scores_sorted,scores),axis=1).T]).T
    
    # equivalent to: ((scores_greater + 1) / (len(scores_greater) + 1))
    pvals = 1 - ((scores_smaller) / (scores.shape[0] + 1))
    
    # Create DFs
    # pvals = pd.DataFrame(pvals, idx_labels, col_labels)
    
    return pvals
