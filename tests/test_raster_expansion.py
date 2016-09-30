import numpy as np
import numpy.ma as ma
import rasterio

value_raster = "analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd.rank.compressed.tif"
mask_raster = "analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es.rank.compressed.tif"
expand_raster = "analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd.rank_es_matched.compressed.tif"

print(" Reading rasters...")
value_raster = rasterio.open(value_raster)
mask_raster = rasterio.open(mask_raster)
expand_raster = rasterio.open(expand_raster)

print(" Extracting values and masks...")
value_raster_src = value_raster.read(1, masked=True)
value_raster_mask = value_raster_src.mask
n_value_raster_mask = np.sum(value_raster_mask)

mask_raster_src = mask_raster.read(1, masked=True)
mask_raster_mask = mask_raster_src.mask
n_mask_raster_mask = np.sum(mask_raster_mask)

expand_raster_src = expand_raster.read(1, masked=True)
expand_raster_mask = expand_raster_src.mask
n_expand_raster_mask = np.sum(expand_raster_mask)

print("\n # of NoData in value_raster mask: \t{}".format(n_value_raster_mask))
print(" # of NoData in mask_raster mask: \t{}".format(n_mask_raster_mask))
print(" # of NoData in expand_raster mask: \t{}".format(n_expand_raster_mask))

print("\n Diff mask_raster and expand_raster: \t{}".format(n_mask_raster_mask - n_expand_raster_mask))

diff = mask_raster_mask != expand_raster_mask
diff_indices = np.where(diff == True)
diff_zipper = zip(diff_indices[0], diff_indices[1])
diff_index_list = [item for item in diff_zipper]

print("")

for coords in diff_index_list:
    x = coords[0]
    y = coords[1]
    vv = value_raster_src[x, y]
    mv = mask_raster_src[x, y]
    ev = expand_raster_src[x, y]
    print( " {0} \t> value_raster: {1:1.10f} | mask_raster: {2} | expand_raster: {3:1.10f}".format(coords, vv, mv, ev))
