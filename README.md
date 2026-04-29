[![PyPI version shields.io](https://img.shields.io/pypi/v/cellex.svg)](https://pypi.python.org/pypi/cellex/)

# CELLEX
CELLEX (CELL-type EXpression-specificity) is a tool for computing cell-type Expression Specificity (ES) profiles. It employs a "wisdom of the crowd"-approach by integrating multiple ES metrics, thus combining complementary cell-type ES profiles, to capture multiple aspects of ES and obtain improved robustness.

![CELLEX_overview](https://user-images.githubusercontent.com/5487016/72679348-9662cf80-3aae-11ea-9d07-c4cea1daec5f.png)

# Contents
- [Documentation](https://github.com/perslab/CELLEX#documentation)
- [Quick Start](https://github.com/perslab/CELLEX#quick-start)
- [Tutorials](https://github.com/perslab/CELLEX#tutorials)
- [Contact and References](https://github.com/perslab/CELLEX#about)

# Documentation
The documentation for CELLEX can be accessed in the following ways:

- **[CELLEX Wiki](https://github.com/perslab/CELLEX/wiki)** : main documentation on the usage of CELLEX
- **[CELLEX API docs](https://perslab.github.io/CELLEX/)**: documentation of CELLEX API/functions
- [**Publication**](https://elifesciences.org/articles/55851): technical details on the CELLEX method. _Genetic mapping of etiologic brain cell types for obesity_ 
([Timshel eLife, 2020](https://elifesciences.org/articles/55851), Appendix)

We are continually updating the documentation for CELLEX. If some information is missing, please submit your request or question via our [issue tracker](https://github.com/perslab/CELLECT/issues).


# Quick start
This brief tutorial showcases the core features of CELLEX.

## TL;DR
```python
import numpy as np
import pandas as pd
import cellex

data = pd.read_csv("./data.csv", index_col=0)
metadata = pd.read_csv("./metadata.csv", index_col=0)

eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)
eso.results["esmu"].to_csv("mydataset.esmu.csv.gz")
```

## Walkthrough
### Setup
#### Install the version from this repo
Clone the development repo and install from source using `pip`. The development version may contain bug fixes that have not been released, as well as experimental features.

```
git clone https://github.com/Tobi1kenobi/CELLEX.git --branch master --single-branch
cd CELLEX
pip install -e .
```

### Import modules
```python
import numpy as np # needed for formatting data for this tutorial
import pandas as pd # needed for formatting data for this tutorial
import cellex
```

### Load input data and metadata
```python
data = pd.read_csv("./data.csv", index_col=0)
metadata = pd.read_csv("./metadata.csv", index_col=0)
```

#### Data format
Data may consist of UMI counts (integer) for each **gene** and **cell**.

|               | cell_1                | ... | cell_9                 |
|---------------|-----------------------|-----|------------------------|
| gene_x        | 0                     | ... | 4                      |
| ...           | ...                   | ... | ...                    |
| gene_z        | 3                     | ... | 1                      |

Shape: *m* genes by *n* cells.

#### Metadata format
Metadata should consist of *unique* cell id's and matching annotation (string).

| cell_id                | cell_type |
|------------------------|-----------|
| cell_1                 | type_A    |
| ...                    | ...       |
| cell_9                 | type_C    |

Shape: *n* cells by 2.

### Create ESObject and compute ESmu

```python
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)

eso.compute(verbose=True)
```

### View Expression Specificity scores
All results are accessible via the `results` attribute of the `ESObject`.

```python
eso.results["esmu"]
```

### Save result(s)
#### Pro-tip: Using CELLEX with CELLECT
The ESmu scores may be used with **[CELLECT](https://github.com/perslab/CELLECT)**. CELLECT requires that genes are in the *Human Ensembl Gene ID* format. CELLEX provides a simple renaming utility for this purpose:

```python
cellex.utils.mapping.mouse_ens_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
```

#### Save ESmu

```python
eso.results["esmu"].to_csv("mydataset.esmu.csv.gz")
```

#### Save all or specific results

```python
eso.save_as_csv(keys=["all"], verbose=True)
```

#### Output format
Output consist of Expression Specificity Weights (float) for each **gene** and **cell-type**. ESmu values lie in the range [0,1].

|               | type_A                | ... | type_C                 |
|---------------|-----------------------|-----|------------------------|
| gene_x        | 0.0                   | ... | 0.9                    |
| ...           | ...                   | ... | ...                    |
| gene_z        | 0.1                   | ... | 0.2                    |

Shape: *m* genes by *x* unique annotations. N.B. a number of genes may be removed during preprocessing.



# Tutorials
Various tutorials and walkthroughs will be made available here, while the Wiki is in the making. These Jupyter Notebooks cover everything from downloading CELLEX and data to analysis and plotting.

* [Demo: Downloading and running CELLEX with sample Mousebrain Atlas data](tutorials/demo_mousebrain_vascular_cells.ipynb)
* [Demo: Downloading and running CELLEX with sample MOCA data](tutorials/demo_moca_100k.ipynb)


# About

## Developers
- Tobias Overlund Stannius (University of Copenhagen) [@TobiasStannius](https://twitter.com/TobiasStannius)
- Pascal Nordgren Timshel (University of Copenhagen) [@ptimshel](https://twitter.com/ptimshel)

## Contact
Please create an [issue](https://github.com/perslab/CELLECT/issues) in this repo, if you encounter any problems using CELLEX. Alternatively, you may write an email to timshel(at)sund.ku.dk

## References

If you find CELLEX useful for your research, please consider citing: 
**[Timshel (eLife, 2020): _Genetic mapping of etiologic brain cell types for obesity_](https://elifesciences.org/articles/55851)**
