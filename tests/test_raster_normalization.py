import numpy as np
import pytest
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


#def test_raster_normalization():

#    assert 11 == 12
