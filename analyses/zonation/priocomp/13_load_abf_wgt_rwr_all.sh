#!/bin/sh
zig4 -l../../RWR/rwr_all_weights_expanded.tif 13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all.dat 13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all.spp 13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all_out/13_load_abf_wgt_rwr_all.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
