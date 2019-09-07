import numpy as np
import pandas as pd
import time
import datetime

def remove_non_expressed(df: pd.DataFrame, verbose: bool=False):
    """
    Identifies and removes genes with mean 0 from the dataframe,
    i.e. non-expressed genes.
    
    Parameters
    ----------
    df : DataFrame
        Gene expression data, e.g. UMI counts.
    
    verbose : bool, optional (default: False)
        Print progress report.
    
    Returns
    -------
    df_filtered : DataFrame
        Filtered data with no genes with mean 0 expression across all cells.

    NOTE:
    * df format – rows: genes, cols: cells.
              A     B     C
        POMC  0.0   0.5   0.9
        AGRP  0.2   0.0   0.0
        LEPR  0.1   0.1   0.4

    TODO:
    * consider a faster way of determining this. Perhaps just check if sum > 0?
    * verbose printing
    """
    start = 0

    if verbose:
        start = time.time()
        print("Preprocessing - running remove_non_expressed ... ", end='')
    
    mask = (df.sum(axis=1) / df.shape[1]) != 0

    df_filtered = df[mask]

    if verbose:
        n_genes_org = len(df)
        n_genes_filtered = len(df_filtered)
        td = datetime.timedelta(seconds=(time.time() - start))
        print("excluded %d / %d genes in %d min %d sec" % ((n_genes_org - n_genes_filtered), n_genes_org, *divmod(td.seconds, 60)))

    return df_filtered
