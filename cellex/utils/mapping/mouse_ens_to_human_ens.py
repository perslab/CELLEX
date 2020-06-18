import numpy as np
import pandas as pd
import pkg_resources

def mouse_ens_to_human_ens(df_unmapped: pd.DataFrame, drop_unmapped: bool=False, verbose: bool=False) -> None:
    """Map dataframe index mouse ensemble id to human ensemble id inplace

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
    
    PREFIX = "ENSMUSG"
    
    if verbose:
        print("Mapping: mouse ensembl gene id's --> human ensembl gene id's ...")
    
    # Check that genes are correct format
    mask_peek = np.array([PREFIX in str(idx) for idx in df_unmapped.index.values])

    if not (mask_peek.any()):
        print("Dataframe index contains values that are not ensemble format or not mouse ensembl id: ", df_unmapped.index.values[mask_peek])
    
    resource_package = __name__
    resource_path = 'maps/hsapiens_mmusculus_unique_orthologs.GRCh38.ens_v100.txt.gz'  # Do not use os.path.join()
    resource_stream = pkg_resources.resource_stream(resource_package, resource_path)    
    df_map = pd.read_csv(resource_stream, compression='gzip', delim_whitespace=True)
    # create dictionary for mapping
    map_dict = dict(zip(df_map["mmusculus_homolog_ensembl_gene"].ravel(), \
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
