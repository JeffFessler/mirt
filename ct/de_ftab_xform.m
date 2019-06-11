  function fs = de_ftab_xform(ftab, fm)
%|function fs = de_ftab_xform(ftab, fm)
%|
%| convert nonlinear BH function "f" to "f*" using jacobian
%| in
%|	ftab	struct	from de_ftab_build()
%|	fm	[n? M]	f_m samples
%| out
%|	fs	[n? L]	f^*_l samples
%|
%| Copyright 2002-2-15, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

[LL MM] = size(ftab.T);
dim = size(fm);
if MM == 1
	dim = [dim 1];
end
if dim(end) ~= MM, keyboard, error 'bad input dim', end
fm = reshape(fm, prod(dim(1:end-1)), MM); % [*n M]
fs = fm * ftab.T'; % [*n L]
fs = reshape(fs, [dim(1:end-1) LL]); % [n? L]
