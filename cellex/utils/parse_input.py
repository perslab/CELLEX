import numpy as np
import pandas as pd
import datetime
import time

def parse_input(data: pd.DataFrame, annotation, verbose: bool=False) -> None:
    """Parse input data and annotation
    
    Checks that data and annotation / metadata match.
    Checks for duplicates.
    Appends _dp{n} to duplicated cell id's.
    
    
    Parameters
    ----------
    data : DataFrame
        Original expression data.
    
    annotation : ndarray, Series, DataFrame
        Annotation to group cells by.
    
    verbose : bool, optional (default: False)
        Print progress report.
    
    Returns
    -------
    data, annotation : (DataFrame, ndarray)
        Data and annotation with checked index/column headers
    
    """
    start = 0

    if verbose:
        start = time.time()
        print("Preprocessing - checking input ... ", end='')
    
    # Check input length
    assert (len(annotation) == data.shape[1]), "Number of annotations do not match number of cells."
    
    # If Series, relevant to check if index matches data cols
    if type(annotation) is pd.Series:
        assert all(np.sort(data.columns.values) == np.sort(annotation.index.values)), "Data columns and annotation index values do not match 1:1."

    # Turn annotation into Series
    # dataframe --> series
    if type(annotation) is pd.DataFrame:
        assert (annotation.shape[1] == 1), "DataFrame annotation had unexpected number of columns {}".format(annotation.shape[1])
        annotation = annotation.iloc[:,0]
    
    # numpy array --> series
    if type(annotation) is np.ndarray:
        annotation = pd.Series(data=annotation, index=data.columns.values)
    
    # Handle duplicates
    if any(data.columns.duplicated()):
        
        if verbose:
            print("\n  duplicate cell id's detected ... ", end='')
        
        cols = pd.Series(data.columns)
        dups = cols[cols.duplicated()].unique()
        n_dups = sum([i in dups for i in cols.values])
        
        # iterate over dups and append suffix _dp{0} .. _dp{n}
        for d in dups:
            # check for differing type-annotations
            types = annotation.values[annotation.index == d]
            if verbose and (len(np.unique(types)) > 1):
                print("\n  duplicated id {uid} has >1 type-annotation: {t}".format(uid=d, t=types), end='')
            
            # add suffix
            mask = cols == d
            names = [d + '_dp' + str(i) for i in range(sum(mask))]
            data.columns.values[mask] = names
            annotation.index.values[mask] = names
            
        if verbose:
            print("\n  {n} duplicate id's renamed ... ".format(n=n_dups), end='')

    # Reduce annotation to sorted ndarray
    # series --> numpy array
    if type(annotation) is pd.Series:
        annotation = data.columns.map(annotation, na_action="ignore").values.astype(str)

    if verbose:
        td = datetime.timedelta(seconds=(time.time() - start))
        print("input parsed in %d min %d sec" % (divmod(td.seconds, 60)))
        
    return data, annotation
