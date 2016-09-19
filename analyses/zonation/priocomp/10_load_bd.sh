#!/bin/sh
zig4 -l06_abf_es/06_abf_es_out/06_abf_es.rank.compressed.tif 10_load_bd/10_load_bd.dat 10_load_bd/10_load_bd.spp 10_load_bd/10_load_bd_out/10_load_bd.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
