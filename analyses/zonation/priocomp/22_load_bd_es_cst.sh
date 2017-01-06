#!/bin/sh
zig4 -l14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank_es_matched.compressed.tif 22_load_bd_es_cst/22_load_bd_es_cst.dat 22_load_bd_es_cst/22_load_bd_es_cst.spp 22_load_bd_es_cst/22_load_bd_es_cst_out/22_load_bd_es_cst.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png
