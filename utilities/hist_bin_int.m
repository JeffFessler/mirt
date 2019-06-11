  function [nk center] = hist_bin_int(data, varargin)
%|function [nk center] = hist_bin_int(data, varargin)
%|
%| fast histogram of multidimensional data into equally-spaced bins
%| todo: time this vs accumarray
%|
%| in
%|	data	[N M]		data values to be binned (M-dimensional)
%|				must be integer valued >= 1
%|
%| option
%|	show	1|0		if no nargout, then display/print (default 1)
%|
%| out
%|	nk	[[ncent]]	histogram values: sum(nk(:)) = N
%|	center	{ncent}		cell array of bin centers for each dimension
%|
%| Copyright 2010-05-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(data, 'test'), hist_bin_int_test, return, end
if nargin < 1, ir_usage, end

arg.show = 1;
arg.how = 'accumarray'; % found to be faster than 'histc' option
arg = vararg_pair(arg, varargin);

M = size(data,2);

ncent = nan(1,M);
center = cell(M,1);

if any(data(:) ~= int32(data(:)))
	fail 'not int'
end
data = int32(data); % 2010-08-13 to avoid overflow with int16

switch arg.how
case 'accumarray'

	for id=1:M
		tmp = data(:,id);
		dmin = min(tmp);
		dmax = max(tmp);
		data(:,id) = tmp - dmin + 1;
		center{id} = dmin:dmax;
		ncent(id) = length(center{id});
	end
	nk = accumarray(data, ones(size(data,1),1));

case 'histc'

	tmp = data(:,1);
	dmin = min(tmp);
	dmax = max(tmp);
	list = tmp - dmin;
	center{1} = dmin:dmax;
	ncent(1) = numel(center{1});
	nprod = ncent(1);

	for id=2:M
		tmp = data(:,id);
		dmin = min(tmp);
		dmax = max(tmp);
		center{id} = dmin:dmax;
		list = list + (tmp - dmin) * nprod;
		ncent(id) = length(center{id});
		nprod = nprod * ncent(id);
	end

	if nprod > 2^30-1, fail 'int32 too small', end

	nk = histc(list, 0:(nprod-1)); % note: specifies left "edges"
	nk = reshape(nk, [ncent 1]);

case 'sparse' % obsolete

	for id=1:M
		tmp = data(:,id);
		dmin = min(tmp);
		dmax = max(tmp);
		data(:,id) = tmp - dmin;
		center{id} = dmin:dmax;
		ncent(id) = length(center{id});
	end

	% the sparse trick below does the following incrementing:
	% for ii=1:N, nk(list(ii)) += 1, end
	s = cumprod(ncent);
	s = [1 s(1:end-1)];
	list = 1 + double(data) * s(:); % need double to avoid overflow
	N = size(data,1);
	s = sparse((1:N)', list, ones(N,1), N, prod(ncent));
	nk = full(sum(s));
	nk = reshape(nk, [ncent 1]);

otherwise
	fail('how=%s not done', arg.how)
end


if ~nargout && arg.show && M <= 2
	if length(nk) < 10
		for ik=1:length(nk)
			printf('hist[%g]=%d', center{1}(ik), nk(ik))
		end
	else
		if M == 1
			plot(center{1}, nk)
		elseif M == 2
			im(center{1}, center{2}, nk)
		end
	end
clear center nk
end


%
% self test
%
function hist_bin_int_test
data = [2 4; 1 2; 2 4];
nk = hist_bin_int(data);
ideal = zeros(2,3);
m1 = min(data(:,1)) - 1;
m2 = min(data(:,2)) - 1;
ideal(2-m1,4-m2) = 2;
ideal(1-m1,2-m2) = 1;
jf_equal(nk, ideal)

data = [3 3 2 5]';
[nk center] = hist_bin_int(data);
jf_equal(nk, [1 2 0 1]')
jf_equal(center{1}, [2 3 4 5])

% look at joint gaussian distribution
rng(0)
n = 2^19;
sig = 5;
rho = 0.7; % correlated gaussian
Cov = sig * [1 rho; rho 1];
tmp = sqrtm(Cov);
data = randn(n, 2) * sqrtm(Cov);
data = round(3 * (data + 8));

types = {'accumarray', 'histc', 'sparse'};
for it=1:numel(types)
	hist_bin_int(data, 'how', types{it}); % warm up

	cpu etic
	[tmp cent] = hist_bin_int(data, 'how', types{it});
	cpu('etoc', types{it})
	if it == 1
		nk = tmp;
	else
		jf_equal(tmp, nk)
	end
end

if im
	im plc 1 2
	im(1, cent{1}, cent{2}, nk)
	axis equal
	im subplot 2
	plot(data(:,1), data(:,2), '.')
	axis equal
end

% check marginals
if 1
	pj = nk/sum(nk(:));
	p1 = sum(pj, 2);
	p2 = sum(pj, 1)';

	n1 = hist_bin_int(data(:,1)) / n;
	n2 = hist_bin_int(data(:,2)) / n;
	equivs(n1, p1)
	equivs(n2, p2)
end
