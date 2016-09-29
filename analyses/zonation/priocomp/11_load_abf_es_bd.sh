#!/bin/sh
zig4 -l06_abf_es/06_abf_es_out/06_abf_es.rank_bd_matched.compressed.tif 11_load_abf_es_bd/11_load_abf_es_bd.dat 11_load_abf_es_bd/11_load_abf_es_bd.spp 11_load_abf_es_bd/11_load_abf_es_bd_out/11_load_abf_es_bd.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
