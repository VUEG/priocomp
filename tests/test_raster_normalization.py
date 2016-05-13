import numpy as np
import numpy.testing as npt
import os
import py
import pytest
import rasterio
import scipy.stats
from src.data_processing.python.rescale import normalize, rescale_raster, standardize

nrow = 100
ncol = 100

# Create random data
np.random.seed(0)

# Floats
positive_float_data_0_1 = np.random.uniform(low=0.0, high=1.0, size=(nrow, ncol))
positive_float_data_0_100 = np.random.uniform(low=0.0, high=100.0, size=(nrow, ncol))
negative_float_data_0_1 = np.random.uniform(low=-1.0, high=0.0, size=(nrow, ncol))
negative_float_data_0_100 = np.random.uniform(low=-100.0, high=0.0, size=(nrow, ncol))

# Integers
positive_int_data_0_100 = np.random.random_integers(low=0, high=100, size=(nrow, ncol))
negpos_int_data_100_100 = np.random.random_integers(low=-100, high=100, size=(nrow, ncol))
negative_int_data_0_100 = np.random.random_integers(low=-100, high=0, size=(nrow, ncol))


def test_normalization_value_types():
        with pytest.raises(TypeError):
            normalize("foo")

def test_normalization_value_range():
    norm_values = normalize(positive_float_data_0_1)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(positive_float_data_0_100)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(negative_float_data_0_1)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(negative_float_data_0_100)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(positive_int_data_0_100)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(negpos_int_data_100_100)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)

    norm_values = normalize(negative_int_data_0_100)
    assert (norm_values.all() >= 0.0) & (norm_values.all() <= 1.0)


def test_normalization_value_rank():
    # Test if the original ranks of values are retained

    def compare_ranks(x, y):
        org_rank = scipy.stats.rankdata(x, method='ordinal')
        new_rank = scipy.stats.rankdata(y, method='ordinal')
        return (org_rank == new_rank)

    assert compare_ranks(positive_float_data_0_1, normalize(positive_float_data_0_1))
    assert compare_ranks(positive_float_data_0_100, normalize(positive_float_data_0_100))
    assert compare_ranks(negative_float_data_0_1, normalize(negative_float_data_0_1))
    assert compare_ranks(negative_float_data_0_100, normalize(negative_float_data_0_100))
    assert compare_ranks(positive_int_data_0_100, normalize(positive_int_data_0_100))
    assert compare_ranks(negpos_int_data_100_100, normalize(negpos_int_data_100_100))
    assert compare_ranks(negative_int_data_0_100, normalize(negative_int_data_0_100))


def test_raster_rescaling_args(tmpdir):
    # Input must be found
    correct_input_filepath = os.path.join(os.path.dirname(__file__), 'data/test_positive_float_data_0_1.tif')
    faulty_input_filepath = os.path.join(os.path.dirname(__file__), 'data/test_positive_float_data_XXX.tif')
    with pytest.raises(OSError):
        rescale_raster(faulty_input_filepath, "data/out.tif", method="normalize")

    with pytest.raises(TypeError):
        rescale_raster(correct_input_filepath, "data/out.tif", method="average")


def test_raster_rescaling_read_write(tmpdir):
    # Test reading in the data from disk, rescaling, writing, and re-reading

    def read_write_reread(x_file):
        with rasterio.open(x_file) as src:
            x = src.read(1)
            # Save raster in a temporary dir
            outfile = tmpdir.mkdir("sub").join(os.path.basename(x_file))
            rescale_raster(x_file, str(outfile), method="normalize")

            # Re-read the raster and see if it's the same
            with rasterio.open(str(outfile)) as new_src:
                new_x = new_src.read(1)
                assert (x == new_x).all()

    read_write_reread(os.path.join(os.path.dirname(__file__), 'data/test_positive_float_data_0_1.tif'))
