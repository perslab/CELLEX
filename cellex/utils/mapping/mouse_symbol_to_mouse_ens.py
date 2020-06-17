import numpy as np
import pandas as pd
import pkg_resources

def mouse_symbol_to_mouse_ens(df_unmapped: pd.DataFrame, drop_unmapped: bool=False, verbose: bool=False) -> None:
    """Map dataframe index mouse gene symbol to mouse ensemble id inplace

    Parameters
    ----------
    df_unmapped : DataFrame
        DataFrame with gene names of appropriate format (e.g. ensembl) as index
    
    drop_unmapped : bool, optional (default: False)
        Remove unmapped genes (rows) from df, False: keep original index.
    
    verbose : bool, optional (default: False)
        Print progress report.

    Returns
    -------
        None
    
    TODO
    ----
        * make one mapping-function for all cases
        * modify drop_unmapped to unmapped: {"drop", "keep", "na"}
        * support for custom mapping file
        * handle lower/upper/mixed casee letters in index
    
    """

    assert (len(df_unmapped) > 0), "Empty dataframe."
   
    if verbose:
        print("Mapping: mouse gene symbols --> mouse ensembl gene id's ...")
    
    resource_package = __name__
    resource_path = 'maps/Mus_musculus.GRCm38.ens_v90.gene_name_version2ensembl.txt.gz'  # Do not use os.path.join()
    resource_stream = pkg_resources.resource_stream(resource_package, resource_path)    
    df_map = pd.read_csv(resource_stream, compression='gzip', delim_whitespace=True)
    # create dictionary for mapping
    map_dict = dict(zip(df_map["gene_name_optimal"].ravel(), \
                            df_map["ensembl_gene_id"].ravel()))
    
    # map genes in-place,
    # i.e. indexes are replaced directly in df
    df_unmapped.rename(index=map_dict, inplace=True)
    
    if verbose or drop_unmapped:
        # check for unmapped genes
        # note the tilde ~ to get genes NOT mapped
        mask_unmapped = ~df_unmapped.index.isin(df_map["ensembl_gene_id"])
        label_unmapped = df_unmapped.index.values[mask_unmapped]
    
        # create report
        n_unmapped = len(label_unmapped)
        
        if verbose:
            n_total = len(df_unmapped)
            pct = n_unmapped / n_total * 100
            print("%.2f pct of genes are unmapped ..." % pct)
        
        if drop_unmapped:
            df_unmapped.drop(index=label_unmapped, inplace=True)
            n_mapped = len(df_unmapped)
            if verbose:
                print("Removed {} unmapped genes ...".format(n_unmapped))
    
    return None
