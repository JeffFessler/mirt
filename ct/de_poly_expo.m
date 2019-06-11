  function expo = de_poly_expo(MM, varargin)
%|function expo = de_poly_expo(MM, [option])
%|
%| Generate exponents for multi-dimensional polynomials with terms:
%|	x[n,1]^expo(k,1) * ... * x[n,M]^expo(k,M)
%|
%| in
%|	MM			dimension of polynomial
%|
%| option
%|	'emax'	[1]|[M]		maximum exponent of any given term. default: 1
%|	'esum'	[1]		maximum of *sum* of exponents. default: inf
%|	'sort'	1|0		sort in terms of increasing power (sum).  default: 1
%|
%| out
%|	expo	[K M]		exponents for using in de_poly_eval()
%|
%| Copyright 2008-10-28, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(MM, 'test'), de_poly_expo_test, return, end

arg.emax = 1;
arg.esum = inf;
arg.sort = true;
arg = vararg_pair(arg, varargin);

if isscalar(arg.emax)
	arg.emax = repmat(arg.emax, [1 MM]);
end

% powers for each term
for mm=1:MM
	expo{mm} = [0:arg.emax(mm)];
end

expo = ndgrid_jf('mat', expo{:});
expo = reshapee(expo, [], MM);
expo = expo(sum(expo,2) <= arg.esum,:);

if arg.sort
	[tmp ii] = sort(sum(expo, 2));
	expo = expo(ii,:);
end


%
% de_poly_expo_test()
%
function de_poly_expo_test

expo = de_poly_expo(3);
expo = de_poly_expo(2, 'emax', 3, 'esum', 4);
