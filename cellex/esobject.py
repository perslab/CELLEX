import h5py
import numpy as np
import os
import pandas as pd
import sys
from .summarydata import SummaryData
from . import preprocessing
from . import metrics
from . import utils
from cellex import ES_METRICS


class ESObject(object):
    """A class that integrates the CELLEX workflow to compute ES

    Workflow: filter non-expressed, normalize, filter non-varying, compute
    ESw, compute p-values, compute ESw*, compute ESmu and ESsd


    Attributes
    ----------
    results : dict
    summary_data : SummaryData

    Methods
    -------
    compute(self, esms: list=None, verbose: bool=False, compute_meta: bool=True) -> None:
    save_as_csv(self, file_prefix: str=None, path: str=None, keys: list=None, verbose: bool=False) -> None:
    save_as_hdf(self, filename: str, path: str=None, keys: list=None, verbose: bool=False) -> None:
    """
    
    def __init__(self, 
                data: pd.DataFrame, 
                annotation: pd.Series,
                dtype: str="float32",
                remove_non_expressed: bool=True,
                normalize: bool=True,
                anova: bool=True,
                verbose: bool=False
                ):
        """
        Parameters
        ----------
        data: DataFrame
            Original expression data.
        
        annotation : Series
            Annotation to group cells by.
        
        dtype : str, optional (default: "float32)
            Datatype to use.
        
        remove_non_expressed : bool, optional (default: True)
            Remove non-expressed genes from dataframe.
        
        normalize : bool, optional (default: True)
            Normalize data.
        
        anova : bool, optional (default: True)
            Perform anova to remove genes with low variance.
        """
        
        ### Preprocessing steps
        # default: cast input datatype to type with smaller memory footprint
        data = data.astype(dtype=dtype)

        # parse data and metadata, i.e. run various quality checks
        data, annotation = utils.parse_input(data, annotation, verbose)
        
        if remove_non_expressed:
            data = preprocessing.remove_non_expressed(df=data, verbose=verbose)
        
        if normalize:
            data = preprocessing.log_normalize(df=data, verbose=verbose)
        
        if anova:
            # anova returns dict of two dataframes. Select the filtered "df"
            data = preprocessing.anova(df=data, annotation=annotation, verbose=verbose)["df"]

        self.results = {}

        ### Create SummaryData object
        self.summary_data = SummaryData(data, annotation)
    
    def compute(self,
                esms: list=None, \
                verbose: bool=False, \
                compute_meta: bool=True) -> None:
        """Compute ESw using specified ES metrics using object annotation.

        Results are stored in a dictionary with key:
            <esm>.<esw_type>

        Parameters
        ----------
        esms : list(str), optional (default: None)
            List of ES metrics to compute
        verbose : bool, optional (default: False)
            Print progress report.
        compute_meta: bool, optional (default: True)
            Compute ESmu, ESsd, and prerequisites.

        Returns
        -------
            None


        """
        
        if esms is None:
            esms = ES_METRICS

        for m in esms:
            if (m.lower() not in ES_METRICS):
                raise ValueError("No such metric: ", m.lower())

        results = {}

        for m in esms:
            esm_result = getattr(metrics, m.lower())(self.summary_data, verbose, compute_meta)
            results.update(esm_result)

        if compute_meta:
            esws = [val for key,val in results.items() if ("esw_s" in key)]
            results["esmu"] = metrics.es_mu(esws, verbose)
            results["essd"] = metrics.es_sd(esws, verbose)

        if verbose:
            print("Computed %a." % list(results.keys()))

        self.results.update(results)

        return None

    def save_as_csv(self, file_prefix: str=None, path: str=None, keys: list=None, verbose: bool=False) -> None:
        """Save results as multiple csv files

        Saves all results in self.results to a directory: /out_###
        Results include esw, esw_null, pvals, qvals and ESmu.

        Parameters
        ----------
        file_prefix : str, optional (default: None)
            Prefix to append to filenames, i.e. 
            <prefix>.<metric>.<item>.csv.gz

        path : str, optional (default: None)
            Path to save to. If None, saves to "out".

        keys : list, optional (default: all keys in self.results)
            Keys of results in self.results dictionary.
            May be used for saving only specific results.

        verbose : bool, optional (default: False)
            Print progress report.

        Returns
        -------
        None
        """
        if verbose:
            print("Saving results as csv to disk ...")
        
        if keys is None:
            keys = ["esmu", "essd"]

        if "all" in keys:
            keys = self.results.keys()

        if path is None:
            path = "out"

        if file_prefix is None:
            file_prefix = ""
        else:
            file_prefix = "{}.".format(file_prefix)

        os.makedirs(path, exist_ok=True) # make dir if it doesn't already exist

        ### Save results
        for k in keys:
            try:
                df = self.results[k]
                fp = "{}/{}{}.csv.gz".format(path, file_prefix, k)
                df.to_csv(fp, compression="gzip")
                if verbose:
                    print("  Saved: {}".format(fp))
            except KeyError:
                print("  WARNING: Key \"{}\" does not exist in ESObject results. No data saved.".format(k))

        if verbose:
            print("Finished saving results to {}".format(path))
    

    def save_as_hdf(self, filename: str, path: str=None, keys: list=None, verbose: bool=False) -> None:
        """Save results to a single hdf file

        Parameters
        ----------
        filename : str
            Filename to write to.

        path : str, optional (default: None)
            Path to save to. If None, saves to "out".

        keys : list, optional (default: all keys in self.results)
            Keys of results in self.results dictionary.
            May be used for saving only specific results.

        verbose : bool, optional (default: False)
            Print progress report.

        Returns
        -------
        None

        """
        if verbose:
            print("Saving results as hdf to disk ...")
        
        if keys is None:
            keys = self.results.keys()
        
        assert (len(keys) > 0), "No results to save in ESObject.results."

        if path is None:
            path = "out"
        
        filename = "{}/{}.h5".format(path, filename)

        os.makedirs(path, exist_ok=True) # make dir if it doesn't already exist

        with h5py.File(filename, "w-") as f:
            axis0 = self.summary_data.mean.columns.values.astype("S")
            axis1 = self.summary_data.mean.index.values.astype("S")
            f["metadata/axis0"] = axis0
            f["metadata/axis1"] = axis1

            if verbose:
                print("  Saved: metadata/axis0")
                print("  Saved: metadata/axis0")

            for k in keys:
                f["data/{}".format(k)] = self.results[k].values
                if verbose:
                    print("  Saved: data/{}".format(k))

            f.flush()

        if verbose:
            print("Finished saving results to {}".format(filename))
