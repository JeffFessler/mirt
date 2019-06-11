 function data = poisson(xm, seed, varargin)
%function data = poisson(xm, seed, [options])
%|
%| option
%|	'factor'	see poisson2.m
%|
%| Generate Poisson random vector with mean xm.
%| For small means, use poisson1.m
%| For large means, use poisson2.m
%| see num. rec. C, P. 222
%|
%| Copyright 1997-4-29, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(xm, 'test'), poisson_test, return, end

if ~isvar('seed')
	seed = [];
end

arg.factor = 0.85;
arg.try_poissrnd = false;
arg = vararg_pair(arg, varargin);

if arg.try_poissrnd % use matlab/octave poissrnd if available
	if exist('poissrnd') == 2
		if ~isempty(seed), rng(seed), end
		data = poissrnd(xm);
		return
	end
end

data = xm;
%xm = xm(:);

small = xm < 12;
data( small) = poisson1(xm( small), seed);
data(~small) = poisson2(xm(~small), 'seed', seed, 'factor', arg.factor);


% poisson_test()
% run a timing race against matlab's poissrnd
function poisson_test
n = 2^8;
t = reshape(linspace(1, 1000, n^2), n, n);
if 1
	poisson(1); % warm-up
	cpu etic
	poisson(t);
	cpu etoc 'fessler poisson time'
end
if exist('poissrnd') == 2
	poissrnd(1); % warm-up
	cpu etic
	poissrnd(t);
	cpu etoc 'matlab poissrnd time'
else
	warn 'matlab poissrnd unavailable'
end

if 0 % look for bad values?
	poisson(10^6 * ones(n,n));
	poisson(1e-6 * ones(n,n));
	poisson(linspace(11,13,2^10+1));
	poisson(linspace(12700,91800,n^3));
end
