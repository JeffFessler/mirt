This directory and ../def/* contain all the source code needed to compile
penalty_mex
for users of unsupported systems, i.e., windows.

Self compiling and self compiled versions are unsupported...

One PC user (thank you!) reported that he was able to get this to compile
from withing matlab using this command:

 mex penalty_mex.c ../def/mexarg.c 'penalty,diff.c' -I../def/ -DMmex -std=c99 -DIs_pc

linux version (should not be needed because i provide mex files for linix):

 mex penalty_mex.c ../def/mexarg.c 'penalty,diff.c' -I../def/ -DMmex CFLAGS='-std=c99 -fPIC'

you will need to move the resulting .mex* file to your matlab path,
e.g., to ../v7

# -outdir ../v7

The "Makefile" is unlikely to be useful.
