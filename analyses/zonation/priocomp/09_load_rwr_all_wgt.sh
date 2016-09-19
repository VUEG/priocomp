#!/bin/sh
zig4 -l../../RWR/rwr_eu26_all_weights.tif 09_load_rwr_all_wgt/09_load_rwr_all_wgt.dat 09_load_rwr_all_wgt/09_load_rwr_all_wgt.spp 09_load_rwr_all_wgt/09_load_rwr_all_wgt_out/09_load_rwr_all_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
