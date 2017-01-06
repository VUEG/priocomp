#!/bin/sh
zig4 -l../../RWR/rwr_all_weights_costs_expanded.tif 25_load_rwr_all_cst/25_load_rwr_all_cst.dat 25_load_rwr_all_cst/25_load_rwr_all_cst.spp 25_load_rwr_all_cst/25_load_rwr_all_cst_out/25_load_rwr_all_cst.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png
