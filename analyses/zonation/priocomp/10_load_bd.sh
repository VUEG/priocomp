#!/bin/sh
zig4 -l08_abf_bd/08_abf_bd_out/08_abf_bd.rank_expanded.compressed.tif 10_load_bd/10_load_bd.dat 10_load_bd/10_load_bd.spp 10_load_bd/10_load_bd_out/10_load_bd.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
