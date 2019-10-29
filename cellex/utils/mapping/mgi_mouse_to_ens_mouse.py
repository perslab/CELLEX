import numpy as np
import pandas as pd
import pkg_resources

def mgi_mouse_to_ens_mouse(df_unmapped: pd.DataFrame, drop_unmapped: bool=False, verbose: bool=False) -> None:
    """
    Maps mouse ensembl gene id's to human ensembl gene id's.

    Args:
        df_unmapped:    a dataframe in tidy-format.
        drop_unmapped:  True: remove unmapped genes (rows) from df, False: keep original index
        verbose:        explicitly print status or not

    Returns:
        None
    
    Todo:
        * modify drop_unmapped to unmapped: {"drop", "keep", "na"}
        * make one mapping-function for all cases
        * support for custom mapping file
        * handle case for empty df
    """

    assert (len(df_unmapped) > 0), "Empty dataframe."
   
    if verbose:
        print("Mapping: mouse mgi gene id's --> mouse ensembl gene id's ...")
    resource_package = __name__
    resource_path = 'maps/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz'  # Do not use os.path.join()
    resource_stream = pkg_resources.resource_stream(resource_package, resource_path)    
    df_map = pd.read_csv(resource_stream, compression='gzip', delim_whitespace=True)
    # create dictionary for mapping mouse ensemble gene id's to human ensembl gene id's
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
