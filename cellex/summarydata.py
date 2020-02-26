import os
import numpy as np
import pandas as pd

class SummaryData(object):
    """A class that contains summary data for computing ES metrics

    The SummaryData object allows sharing of summary data in between 
    commputations, thus reducing the number computations needed,
    e.g. we need only compute the mean once.

    Attributes
    ----------
    data: DataFrame
        Original expression data.

    annotation : 
        Annotation to group cells by.

    _mean : DataFrame
        Mean expression for the groups specified by the annotation.

    _n_cells_per_anno : DataFrame
        Number of cells per groups specified by the annotation.

    _n_nonzero : DataFrame
        Number of nonzero expression values per groups specified by the annotation.

    _variance : DataFrame
        Variance of expression values per groups specified by the annotation.

    annotation_null : np.ndarray
        Null annotation to group cells by.
    
    _mean_null : DataFrame
        Mean expression for the groups specified by the null annotation.

    _n_cells_per_anno_null : DataFrame
        Number of cells per groups specified by the null annotation.

    _n_nonzero_null : DataFrame
        Number of nonzero expression values per groups specified by the null annotation.

    _variance_null : DataFrame
        Variance of expression values per groups specified by the null annotation.

    Methods
    -------
    mean(self)
        Lazily evaluate and return _mean.
    
    n_cells_per_anno(self)
        Lazily evaluate and return _n_cells_per_anno.
    
    n_nonzero(self)
        Lazily evaluate and return _n_nonzero.
    
    variance(self)
        Lazily evaluate and return _variance.
    
    annotation_null(self)
        Lazily evaluate and return _annotation_null.
    
    mean_null(self)
        Lazily evaluate and return _mean_null.

    n_cells_per_anno_null(self)
        Lazily evaluate and return _n_cells_per_anno_null.
    
    n_nonzero_null(self)
        Lazily evaluate and return _n_nonzero_null.
    
    variance_null(self)
        Lazily evaluate and return _variance_null.

    save(dir_name: str=None, verbose: bool=False)
        Save the object attributes, excluding the original data.
    """

    def __init__(self, data: pd.DataFrame, annotation: np.array):
        """
        Parameters
        ----------
        data: DataFrame
            Original expression data.

        annotation : array
            Annotation to group cells by.
        """

        # df.columns = pd.MultiIndex.from_arrays([df.columns,
        #         df.columns.map(annotation, na_action="ignore").values.astype(str)],
        #         names=("id", "annotation"))
        
        self.data = data

        self.annotation = annotation
        self._mean = None
        self._n_cells_per_anno = None # tstat, ges (just one line)
        self._n_nonzero = None
        self._variance = None # tstat
        
        # null attributes
        self._annotation_null = None
        self._mean_null = None
        self._n_nonzero_null = None
        self._n_cells_per_anno_null = None
        self._variance_null = None

    @property
    def mean(self):
        """Compute or return mean"""
        if self._mean is None:
            self._mean = self.data.groupby(self.annotation, axis=1).mean()
        
        return self._mean

    @property
    def n_cells_per_anno(self):
        if self._n_cells_per_anno is None:
            self._n_cells_per_anno = self.data.groupby(self.annotation, axis=1).size()

        return self._n_cells_per_anno

    @property
    def n_nonzero(self):
        if self._n_nonzero is None:
            self._n_nonzero = self.data.groupby(self.annotation, axis=1).agg(lambda x: x.ne(0).sum(axis=1))
        
        return self._n_nonzero

    @property
    def variance(self):
        if self._variance is None:
            self._variance = self.data.groupby(self.annotation, axis=1).var()

        return self._variance

    ### Null summary data
    @property
    def annotation_null(self):
        if self._annotation_null is None:
            np.random.seed(1)
            self._annotation_null = np.random.permutation(self.annotation)

        return self._annotation_null

    @property
    def mean_null(self):
        if self._mean_null is None:
            self._mean_null = self.data.groupby(self.annotation_null, axis=1).mean()

        return self._mean_null
    
    @property
    def n_cells_per_anno_null(self):
        if self._n_cells_per_anno_null is None:
            self._n_cells_per_anno_null = self.data.groupby(self.annotation_null, axis=1).size()

        return self._n_cells_per_anno_null

    @property
    def n_nonzero_null(self):
        if self._n_nonzero_null is None:
            self._n_nonzero_null = self.data.groupby(self.annotation_null, axis=1).agg(lambda x: x.ne(0).sum(axis=1))
        
        return self._n_nonzero_null

    @property
    def variance_null(self):
        if self._variance_null is None:
            self._variance_null = self.data.groupby(self.annotation_null, axis=1).var()

        return self._variance_null

    def save(self, dir_name: str=None, verbose: bool=False) -> None:
        """
        Save summary statistics to disk.
        """
        if verbose:
            print("Saving results to disk ...")

        if dir_name == None:
            dir_name = "out"

        os.makedirs(dir_name, exist_ok=True) # make dir if it doesn't already exist

        ### Loop over SummaryStats attributes and save stats to disk
        for s in dir(self):
            att = getattr(self, s)
            if isinstance(att, pd.DataFrame) and s != "data":
                fp = "{}/summarystat.{}.csv.gz".format(dir_name, self.name, s)
                att.to_csv(fp, compression="gzip")
                if verbose:
                    print("  Saved: {}".format(fp))

            