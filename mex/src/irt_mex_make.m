% irt_mex_make.m

if ispc
	cd mex/src/penalty
    mex penalty_mex.c mexarg.c 'penalty,diff.c' -DIs_pc -DMmex -outdir ../../local
	cd ../../..
else
	cd mex/src/penalty
    mex -v penalty_mex.c '../def/mexarg.c' 'penalty,diff.c' -I'../def/' -DMmex -outdir ../../local
	cd ../../..
end
