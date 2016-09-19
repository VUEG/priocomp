#!/bin/sh
zig4 -l../../ILP/ilp_eu26_all_weights.tif 10_load_ilp_all_wgt/10_load_ilp_all_wgt.dat 10_load_ilp_all_wgt/10_load_ilp_all_wgt.spp 10_load_ilp_all_wgt/10_load_ilp_all_wgt_out/10_load_ilp_all_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
