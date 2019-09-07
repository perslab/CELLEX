import numpy as np
import pandas as pd
import time
import datetime

def es_sd(esws: list, verbose: bool=False) -> pd.DataFrame:
    """Compute ESsd
    
    Computes standard deviation of ESWs, i.e. ES_sd, for each gene / cell-type.

    Parameters
    ----------
    esws : list
        List of ESw dataframes to compute mean from.

    verbose : bool, optional (default: False)
        Print progress report.
    
    Returns
    -------
    es_sd : DataFrame
        ESw std. dev.
    """

    start = 0

    if verbose:
        start = time.time()
        print("Computing ESsd ...")

    if len(esws) < 4:
        print("WARNING: Computing esw_sd using {} metrics ...".format(len(esws)))

    es_sd = pd.DataFrame(data=np.std(([df.values for df in esws]), axis=0), 
                            columns=esws[0].columns.values, 
                            index=esws[0].index.values)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))
    
    return es_sd
