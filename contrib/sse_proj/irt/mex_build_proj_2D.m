% $Id: mex_build_proj_2D.m,v 1.2 2007/08/07 22:02:08 jmh Exp $

dir_current = pwd;
dir_proj_2d = path_find_dir('irt');
dir_proj_2d_src = [dir_proj_2d filesep 'mex' filesep 'src'];
cd(dir_proj_2d_src)

mex -I../include jmh_sse_mex.cpp mexarg.c tomo_mex.cpp tomo_lut.cpp ...
    tomo_strip.cpp mm_malloc.cpp -DMmex -DIs_pc -outdir ../v7

cd(dir_current)
