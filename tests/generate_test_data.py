#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import rasterio

nrow = 100
ncol = 100

# Create random data
np.random.seed(0)

# Floats
positive_float_data_0_1 = np.random.uniform(low=0.0, high=1.0, size=(nrow, ncol))
positive_float_data_0_100 = np.random.uniform(low=0.0, high=100.0, size=(nrow, ncol))

# Integers
positive_int_data_0_100 = np.random.random_integers(low=0, high=100, size=(nrow, ncol))
negpos_int_data_100_100 = np.random.random_integers(low=-100, high=100, size=(nrow, ncol))

# Write the product as a raster band to a new 8-bit file. For
# the new file's profile, we start with the meta attributes of
# the source file, but then change the band count to 1, set the
# dtype to uint8, and specify LZW compression.
profile = rasterio.default_gtiff_profile()
profile.update(dtype=rasterio.float32)

# Write out the data
with rasterio.open('data/test_positive_float_data_0_1.tif', 'w', height=nrow, width=ncol, count=1, **profile) as dst:
    dst.write(positive_float_data_0_1.astype(rasterio.float32), 1)

with rasterio.open('data/test_positive_float_data_0_100.tif', 'w', height=nrow, width=ncol, count=1, **profile) as dst:
    dst.write(positive_float_data_0_100.astype(rasterio.float32), 1)

profile.update(dtype=rasterio.uint8)
with rasterio.open('data/test_positive_int_data_0_100.tif', 'w', height=nrow, width=ncol, count=1, **profile) as dst:
    dst.write(positive_int_data_0_100.astype(rasterio.uint8), 1)

profile.update(dtype=rasterio.int16)
with rasterio.open('data/test_negpos_int_data_100_100.tif', 'w', height=nrow, width=ncol, count=1, **profile) as dst:
    dst.write(negpos_int_data_100_100.astype(rasterio.int16), 1)
