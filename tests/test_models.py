import csv
import pathlib

from xopen import xopen

import gwascatalog.gene
import gwascatalog.gene.models as models
import pytest


def test_version():
    # test importing the package and getting the __version__ variable
    assert getattr(gwascatalog.gene, "__version__", None) is not None


def get_test_data_files():
    parent_dir = pathlib.Path(__file__).parent
    data_dir = parent_dir / "data" / "34662886"  # Backman paper
    sumstat_files = list(data_dir.glob("*.tsv.gz"))

    if not sumstat_files:
        msg = f"""
        No test data files found in {data_dir.resolve()}
        Did you remember to download the test data before starting the test suite? 
        Try running nox -s create_test_data 
        """
        raise RuntimeError(msg)

    return sumstat_files


@pytest.mark.parametrize("file_path", get_test_data_files())
def test_sumstat(file_path):
    with xopen(file_path) as fh:
        genes = list(csv.DictReader(fh, delimiter="\t"))
        for gene in genes:
            # pydantic will raise validation errors
            _ = models.GeneModel(**gene)

