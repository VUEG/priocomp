#!/bin/sh
zig4 -l10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank_bd_matched.compressed.tif 21_load_es_bd_cst/21_load_es_bd_cst.dat 21_load_es_bd_cst/21_load_es_bd_cst.spp 21_load_es_bd_cst/21_load_es_bd_cst_out/21_load_es_bd_cst.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png
