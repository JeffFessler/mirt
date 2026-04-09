% irt_mex_make.m

if ispc
	cd penalty
	mex penalty_mex.c mexarg.c 'penalty,diff.c' -DIs_pc -DMmex -outdir ../v7
	cd ..
else
	cd penalty
	mex -v penalty_mex.c mexarg.c 'penalty,diff.c' -DMmex -outdir ../v7
	cd ..
end
