  function [nk, center] = hist_equal(data, ncent, varargin)
%|function [nk, center] = hist_equal(data, ncent, varargin)
%|
%| fast histogram of multidimensional data into equally-spaced bins
%| todo: accumarray?
%|
%| in
%|	data	[N M]		data values to be binned (M-dimensional)
%|	ncent	[1 M]		# of centroids for each dimension
%| option
%|	'ifsame' char		what to do if all data same along some dimension
%|				'orig' use original ncent values (default)
%|				'1bin' ignore ncent value and use 1 bin
%| out
%|	nk	[[ncent]]	histogram values: sum(nk(:)) = N
%|	center	{ncent}		cell array of bin centers for each dimension
%|
%| Copyright 2004-7-5, Jeff Fessler, University of Michigan

if nargin == 1 && streq(data, 'test'), hist_equal_test, return, end
if nargin < 2, ir_usage, end

arg.ifsame = 'orig'; % default
arg.fudge = 1.001;
arg.dmin = [];
arg.dmax = [];
arg = vararg_pair(arg, varargin);

[N M] = size(data);
if M ~= length(ncent), error 'bad dimensions', end

list = zeros(N,M);
for id=1:M
	nc = ncent(id);
	if ~isempty(arg.dmin)
		dmin = arg.dmin(id);
	else
		dmin = min(data(:,id));
	end
	if ~isempty(arg.dmax)
		dmax = arg.dmax(id);
	else
		dmax = max(data(:,id));
	end
	if dmin == dmax
		switch arg.ifsame
		case 'orig'
			center{id} = dmin + [0:(nc-1)]';
			list(:,id) = 1;
		case '1bin'
			ncent(id) = 1;
			center{id} = dmin;
			list(:,id) = 1;
		otherwise
			fail('option ifzero="%s" unknown', arg.ifzero)
		end
	else
		dmin = dmin * arg.fudge;
		dmax = dmax / arg.fudge;
		center{id} = col(linspace(dmin, dmax, nc));
		ddif = center{id}(2) - center{id}(1);

		ii = 1 + floor((data(:,id) - dmin) / ddif);
		ii = min(ii, nc);
		ii = max(ii, 1);

		list(:,id) = ii;
	end
end

% the sparse trick below does the following incrementing:
% for ii=1:N, nk(list(ii)) += 1, end
s = cumprod(ncent);
s = [1 s(1:end-1)];
list = 1 + (list-1) * s(:);
s = sparse([1:N]', list, ones(N,1), N, prod(ncent));
nk = full(sum(s));
nk = reshape(nk, ncent);


%
% self test
%
function hist_equal_test
rng(0)
r1 = 10 * randn(10^3,1);
r2 = 40 + 10 * randn(10^3,1);
i1 = 5 * randn(10^3,1);
i2 = 7 + 3 * randn(10^3,1);
x = [[r1; r2], 1+[i1; i2]];
[nk cents] = hist_equal(x, [18 15], 'ifsame', '1bin');
if im
	clf, subplot(121), plot(x(:,1), x(:,2), '.')
	c1 = cents{1};
	c2 = cents{2};
%	axis([minmax(creal); minmax(cimag)])
	im(122, c1, c2, nk), axis normal
	[c1 c2] = ndgrid(c1, c2);
	subplot(121)
	hold on
	plot(c1, c2, 'g.')
	hold off
	xlabel 'real', ylabel 'imag'
end
