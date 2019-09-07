import numpy as np
import pandas as pd
import time
import datetime

def log_normalize(df: pd.DataFrame, constant: int=1, scale: int=1e4, verbose: bool=False):
    """
    Normalizes data by scaling and log-transformation
    Default scale is 10,000
    
    Parameters
    ----------
    df : DataFrame
        Raw expression values, e.g. UMI counts.
    
    constant : int, optional (default: 1)
        Constant added before taking log.
    
    scale : float, optional (default: 1e4)
        Scale factor for common transcript count.
    
    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
    df_lognorm : DataFrame
        Log normalized expression values

    NOTE:
    * df format – rows: genes, cols: cells.
              A     B     C
        POMC  0     2     1
        AGRP  3     0     0
        LEPR  1     1     4
    
    TODO:
    * verbose printing
    """
    start = 0

    if verbose:
        start = time.time()
        print("Preprocessing - normalizing data ... ", end='')

    # divide df by column sums to normalize values per cell
    # scale df values by scale (default: 1e4)
    # add 1 to all values in df to avoid div/0 errors in log
    # take log of values in df
    
    df_lognorm = np.log((df / df.sum(axis=0) * scale) + constant)
    
    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("data normalized in %d min %d sec" % (divmod(td.seconds, 60)))

    return df_lognorm
