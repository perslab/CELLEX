import numpy as np
import pandas as pd
import time
import datetime

def es_mu(esws: list, verbose: bool=False) -> pd.DataFrame:
    """Compute ESmu
    
    Computes mean of ESWs, i.e. ES_mu, for each gene / cell-type.

    Parameters
    ----------
    esws : list
        List of ESw dataframes to compute mean from.

    verbose : bool, optional (default: False)
        Print progress report.
    
    Returns
    -------
    esmu : DataFrame
        ESw mean.
    """

    start = 0

    if verbose:
        start = time.time()
        print("Computing ESmu ...")

    if len(esws) < 4:
        print("WARNING: Computing esw_mu using {} metrics ...".format(len(esws)))
    
    esmu = pd.DataFrame(data=np.mean(([df.values for df in esws]), axis=0), 
                            columns=esws[0].columns.values, 
                            index=esws[0].index.values)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("    finished in %d min %d sec" % (divmod(td.seconds, 60)))

    return esmu
