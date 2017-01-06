#!/bin/sh
zig4 -lanalyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank_expanded.compressed.tif 18_load_bd_all_cst/18_load_bd_all_cst.dat 18_load_bd_all_cst/18_load_bd_all_cst.spp 18_load_bd_all_cst/18_load_bd_all_cst_out/18_load_bd_all_cst.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png
