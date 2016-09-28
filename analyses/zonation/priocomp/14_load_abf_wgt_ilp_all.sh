#!/bin/sh
zig4 -l../../ILP/ilp_eu26_all_weights.tif 14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all.dat 14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all.spp 14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all_out/14_load_abf_wgt_ilp_all.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
