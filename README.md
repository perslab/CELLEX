# CELLEX
CELLEX (CELL-type EXpression-specificity) is a tool for computing cell-type Expression Specificity (ES) profiles. It employs a "wisdom of the crowd"-approach by integrating multiple ES metrics, thus combining complementary cell-type ES profiles, to capture multiple aspects of ES and obtain improved robustness.



## Contents
* [Quick start](#Quick-start)
* [Documentation](docs/)


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
import cellex as cellex # needed when importing directly from this repo
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
