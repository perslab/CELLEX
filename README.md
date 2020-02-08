# CELLEX
CELLEX (CELL-type EXpression-specificity) is a tool for computing cell-type Expression Specificity (ES) profiles. It employs a "wisdom of the crowd"-approach by integrating multiple ES metrics, thus combining complementary cell-type ES profiles, to capture multiple aspects of ES and obtain improved robustness.

![CELLEX_overview](https://user-images.githubusercontent.com/5487016/72679348-9662cf80-3aae-11ea-9d07-c4cea1daec5f.png)


See [Timshel (bioRxiv, 2020): Mapping heritability of obesity by cell types](https://www.biorxiv.org/content/10.1101/2020.01.27.920033v1) for further details on the CELLEX method. See also the [CELLEX Wiki](https://github.com/perslab/CELLEX/wiki).

# Quick start
This brief tutorial showcases the core features of CELLEX.

## Setup
### Option A: Install the latest release from PyPi
```
pip install cellex
```

### Option B: Install the development version from this repo
Clone this repo
```
git clone https://github.com/perslab/CELLEX.git --branch develop --single-branch
```
and install it using `pip`
```
cd CELLEX
pip install -e .
```

## Import modules
```python
import numpy as np # needed for formatting data for this tutorial
import pandas as pd # needed for formatting data for this tutorial
import cellex
```

## Load input data and metadata
```python
data = pd.read_csv("./data.csv", index_col=0)
metadata = pd.read_csv("./metadata.csv", index_col=0)
```

### Data format
Data may consist of UMI counts (integer) for each **gene** and **cell**.

|               | cell_1                | ... | cell_9                 |
|---------------|-----------------------|-----|------------------------|
| gene_x        | 0                     | ... | 4                      |
| ...           | ...                   | ... | ...                    |
| gene_z        | 3                     | ... | 1                      |

Shape: *m* genes by *n* cells.

### Metadata format
Metadata should consist of *unique* cell id's and matching annotation (string).

| cell_id                | cell_type |
|------------------------|-----------|
| cell_1                 | type_A    |
| ...                    | ...       |
| cell_9                 | type_C    |

Shape: *n* cells by 2.

## Create ESObject and compute ESmu

```python
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)

eso.compute(verbose=True)
```

## Save result(s)
Only saves ESmu by default. The ESmu specificity scores may be used directly with **[CELLECT](https://github.com/perslab/CELLECT)**.

```python
eso.save_as_csv(verbose=True)
```

### Output format
Output consist of Expression Specificity Weights (float) for each **gene** and **cell-type**. ESmu values lie in the range [0,1].

|               | type_A                | ... | type_C                 |
|---------------|-----------------------|-----|------------------------|
| gene_x        | 0.0                   | ... | 0.9                    |
| ...           | ...                   | ... | ...                    |
| gene_z        | 0.1                   | ... | 0.2                    |

Shape: *m* genes by *x* unique annotations. N.B. a number of genes may be removed during preprocessing.



# Tutorials
Various tutorials and walkthroughs will be made available here, while the Wiki is in the making. These Jupyter Notebooks cover everything from downloading CELLEX and data to analysis and plotting.

* [Demo: Downloading and running CELLEX with sample Mousebrain Atlas data](https://nbviewer.jupyter.org/github/perslab/CELLEX/blob/master/docs/tutorials/demo_mousebrain_vascular_cells.ipynb)
* [Demo: Downloading and running CELLEX with sample MOCA data](https://nbviewer.jupyter.org/github/perslab/CELLEX/blob/master/docs/tutorials/demo_moca_100k.ipynb)


# About

## Developers
- Tobias Overlund Stannius (University of Copenhagen)
- Pascal Nordgren Timshel (University of Copenhagen)

## Contact
Please create an issue in this repo, if you encounter any problems using CELLEX. Alternatively, you may write an email to timshel(at)sund.ku.dk

## References

If you find CELLEX useful for your research, please consider citing: 
**[Timshel (bioRxiv, 2020): _Mapping heritability of obesity by cell types_](https://www.biorxiv.org/content/10.1101/2020.01.27.920033v1)**


