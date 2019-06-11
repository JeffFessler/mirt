 function mse = mri_exp_errs(arg, Eh, Llist, ti, zmap)
%function mse = mri_exp_errs(arg, Eh, Llist, ti, zmap)
% evaluate errors of MRI exponential approximations for a range of L values
% see mri_exp_approx()
% in
%	arg	cell	to be passed to mri_exp_approx
%	Eh	[M,K]	K=N or K<<N
%	Llist	[L,1]	
%	ti	[M,1]	time samples
%	zmap	[N,1]	rate map
% out
%	mse	[L,K]
%
% Copyright 2004-7-5, Jeff Fessler, The University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

nL = length(Llist);
mse = zeros(nL, ncol(Eh));
for ll=1:nL
	L = Llist(ll);
	ticker(arg{1}, ll, nL)
	if size(Eh, 2) == length(zmap) % this one pigs memory
		[B C] = mri_exp_approx(ti, zmap, L, arg{:});
	else % histogram-sized C
		[B C] = mri_exp_approx(ti, zmap, L, 'ctest', 1, arg{:});
	end
	mse(ll,:) = mean(abs(Eh - B * C).^2);
end
