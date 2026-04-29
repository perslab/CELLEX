import pandas as pd
import cellex
from pathlib import Path

# To run tests:
# * activate your venv and ensure cellex is not installed
# * navigate to CELLEX/
# * run `python -m pytest`

BENCHMARK_DIR = Path(__file__).parent / Path('_benchmark/')


def test_recreate_benchmark():
    path_data = str(Path(BENCHMARK_DIR, "zeisel2018_mousebrain_subset_data.csv.gz"))
    path_metadata = str(Path(BENCHMARK_DIR, "zeisel2018_mousebrain_subset_metadata.csv.gz"))
    
    data = pd.read_csv(path_data, index_col=0, compression="gzip")
    metadata = pd.read_csv(path_metadata, index_col=0, compression="gzip")

    eso = cellex.ESObject(data, annotation=metadata, dtype="float64")
    eso.compute()

    for key in eso.results.keys():
        path = str(Path(BENCHMARK_DIR, "mousebrain_subset.{}.csv.gz".format(key)))
        
        benchmark = pd.read_csv(path, index_col=0, compression="gzip")
        result = eso.results[key]
        # Check that results are equal with 14 decimal points of precision
        # i.e. (0.048740023444475165 == 0.04874002344447517) is considered equal
        pd.testing.assert_frame_equal(benchmark, result, atol=1e-14)
