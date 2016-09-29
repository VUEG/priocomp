#!/bin/sh
zig4 -l08_abf_bd/08_abf_bd_out/08_abf_bd.rank_es_matched.compressed.tif 12_load_abf_bd_es/12_load_abf_bd_es.dat 12_load_abf_bd_es/12_load_abf_bd_es.spp 12_load_abf_bd_es/12_load_abf_bd_es_out/12_load_abf_bd_es.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
