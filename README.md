# CELLEX
CELLEX (CELL-type EXpression-specificity) is a tool for computing cell-type Expression Specificity (ES) profiles. It employs a "wisdom of the crowd"-approach by integrating multiple ES metrics, thus combining complementary cell-type ES profiles, to capture multiple aspects of ES and obtain improved robustness.

![CELLEX_overview](https://user-images.githubusercontent.com/5487016/72679348-9662cf80-3aae-11ea-9d07-c4cea1daec5f.png)


# Tutorials
See [tutorials](https://github.com/perslab/CELLEX/tree/master/docs/tutorials) for Jupyter Notebook demos and tutorials.

# Documentation

See ["Mapping heritability of obesity by brain cell types" (Timshel, bioRxiv 2020)](UPDATE_LINK_LATER) for details in the CELLEX method. See also the [CELLEX Wiki](https://github.com/perslab/CELLEX/wiki) for additional documentation.

# Quick start
This brief tutorial showcases the core features of CELLEX.

## Setup

Install CELLEX via pip:
```
pip install cellex
```

<!-- Download this repository and place it in the same directory as the script or Jupyter Notebook you wish to use CELLEX with. -->

## Import modules
```python
import numpy as np
import pandas as pd
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
Metadata should consist of cell id's and matching annotation (string).

| cell_id                | cell_type |
|------------------------|-----------|
| cell_1                 | type_A    |
| ...                    | ...       |
| cell_9                 | type_C    |

Shape: *n* cells by 2.

## Create ESObject and compute ESmu

```python
eso = cellex.ESObject(df=data, annotation=metadata, verbose=True)

eso.compute(verbose=True)
```

## Save result(s)
Only saves ESmu by default. The ESmu specificity scores may be used directly with **[CELLECT](https://github.com/perslab/CELLECT)**.

```python
eso.save(verbose=True)
```

### Output format
Output consist of Expression Specificity Weights (float) for each **gene** and **cell-type**. ESmu values lie in the range [0,1].

|               | type_A                | ... | type_C                 |
|---------------|-----------------------|-----|------------------------|
| gene_x        | 0.0                   | ... | 0.9                    |
| ...           | ...                   | ... | ...                    |
| gene_z        | 0.1                   | ... | 0.2                    |

Shape: *m* genes by *x* unique annotations. N.B. a number of genes may be removed during preprocessing.
