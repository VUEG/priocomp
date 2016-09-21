#!/bin/sh
zig4 -l06_abf_es/06_abf_es_out/06_abf_es.rank_expanded.compressed.tif 09_load_es/09_load_es.dat 09_load_es/09_load_es.spp 09_load_es/09_load_es_out/09_load_es.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
