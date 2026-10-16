% irt_mex_make.m

if ispc
	cd penalty
	mex penalty_mex.c ../def/mexarg.c 'penalty,diff.c' -DIs_pc -DMmex -outdir ../v7
	cd ..
else
	cd penalty
    mex -v penalty_mex.c ../def/mexarg.c 'penalty,diff.c' -I../def/ -DMmex CFLAGS='-std=c99 -fPIC' -outdir ../v7
	cd ..
end
