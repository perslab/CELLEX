import numpy as np
import pandas as pd
import time
import datetime
from scipy import stats

def anova(df: pd.DataFrame, annotation: np.ndarray, threshold: float=0.00001, verbose: bool=False):
    """Apply ANOVA and filter data
    Applies ANOVA to identify lowly / sporadically expressed genes and remove them.
    Values of each gene (row) is split into groups according to the annotation.
    ANOVA is then computed for these groups to get a p-value for this gene.
    The p-value is used to filter out genes with insignificant differences in variance.
    
    Parameters
    ----------
    df : DataFrame
        Log normalized expression data.
    
    annotation : ndarray
        Annotation of cells, e.g. cell-type or class.
    
    threshold : float, optional (default: 0.00001)
        P-value cutoff.

    verbose : bool, optional (default: False)
        Print progress report.
    
    Returns
    -------
    dict
        "anova" : df_anova    : DataFrame of F-values and p-values.
        "df"    : df_filtered : DataFrame of genes with pval < threshold.

    NOTE:
    * df format – rows: genes, cols: cells.
              A     B     C
        POMC  0.0   0.5   0.9
        AGRP  0.2   0.0   0.0
        LEPR  0.1   0.1   0.4

    * annotation format – 1 by n, i.e. one row.
        ["ABC", "ABC", "ABC", "ACBG", "ACBG", ..., "ACMB"]

    TODO:
    * refine grouping. Using np.unique seems unsafe and complex
    * consider using generator instead of list comprehension to save memory and time
    * consider alternative masking/filtering. Is numpy the fastest solution?
    """
    start = 0

    if verbose:
        start = time.time()
        print("Preprocessing - running ANOVA ... ", end='')

    ### Split annotation into groups
    # we use pandas.unique, as it is faster than numpy for big arrays
    # furthermore, it does not sort the groups and thus the original
    # order of indexes is maintained, improving locality
    idx = [np.where(annotation == i)[0] for i in pd.unique(annotation)]
    
    ### Iterate over genes and compute ANOVA for each gene
    def _anova_(row, groups):
        """Apply ANOVA to a vector
        Split the annotation by 
        """
        fval, pval = stats.f_oneway(*[row[g] for g in groups])
        return np.array([fval, pval])

    ### Create dataframe for ANOVA results
    f = np.apply_along_axis(func1d=_anova_, axis=1, arr=df.values, groups=(idx))

    df_anova = pd.DataFrame({"statistic": f[:,0], 
                         "pvalue": f[:,1]}, 
                        index=df.index)

    ### Filter original dataframe
    # Get row-index of genes where pvalue < threshold
    idx_thresh = (df_anova.values[:,1] < threshold)
    df_filtered = df.iloc[idx_thresh]

    if verbose:
        n_genes_org = len(df)
        n_genes_filtered = len(df_filtered)
        td = datetime.timedelta(seconds=(time.time() - start))
        print("excluded %d / %d genes in %d min %d sec" % ((n_genes_org - n_genes_filtered), n_genes_org, *divmod(td.seconds, 60)))
    
    return {'anova': df_anova, 'df': df_filtered}
