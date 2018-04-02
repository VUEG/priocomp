#!/bin/sh
zig4 -l12_abf_bd/12_abf_bd_out/12_abf_bd.rank_es_matched.compressed.tif 20_load_bd_es/20_load_bd_es.dat 20_load_bd_es/20_load_bd_es.spp 20_load_bd_es/20_load_bd_es_out/20_load_bd_es.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png
