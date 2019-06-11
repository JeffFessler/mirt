 function centers = lloyd_max_hist(data, centers, MM, tol, max_iter, chat)
%function centers = lloyd_max_hist(data, centers, MM, tol, max_iter, chat)
%|
%| "improved" version of the lloyd-max algorithm for scalar quantizer design
%| that saves computation by histogramming the data first.
%| Before using this, try using highrate_centers() first!
%|
%| in
%|	data	[N 1]	training data
%|	centers	[K 1]	initial guess of centroids (codebook)
%|	MM		# of histogram bins: K < M < N
%| out
%|	centers [K 1]	final centroids
%|
%| Copyright 2004-7-1, Jeff Fessler, University of Michigan

if nargin == 1 && streq(data, 'test'), lloyd_max_hist_test, return, end
if nargin < 2, ir_usage, end
if nargin < 3, MM = 0; end
if nargin < 4, tol = 1e-3; end
if nargin < 5, max_iter = 40; end
if nargin < 6, chat = 0; end

data = data(:);
Norig = length(data);
if MM
	if ~isreal(data)
		if length(MM) == 1
			MM = [MM MM];
		end
		[wt data] = hist_equal([real(data) imag(data)], MM);
		[dr di] = ndgrid(data(:,1), data(:,2));
		data = dr + 1i * di;
	else
		[wt data] = hist(data, MM);
	end
	data = data(:);
	wt = wt(:);
	% eliminate data bins with 0 wt
	data = data(wt ~= 0);
	wt = wt(wt ~= 0);
else
	wt = ones(size(data));
end

K = length(centers);
if length(unique(centers)) ~= K
	error 'initial centers are not unique'
end
redundant = 0;
if K > length(data)
	if K > Norig
		warning(sprintf('#centroids %d > #data %d!?', K, Norig))
	else
		printf('Warn: #centroids %d > #unique(data) %d.', K, length(data))
	end
	redundant = 1;
end

tol = tol * max(abs(data));
iter = 1;
change = inf;
while iter <= max_iter && change > tol
	if ~redundant && (length(unique(centers)) ~= K)
		warning 'centers became not unique'
	end
	index = quant1_index(data, centers);
	old = centers;
	for kk=1:K
		ik = index == kk;
		if sum(ik)
			if ~sum(wt(ik))
				warning 'bug?'
				keyboard
			end
			centers(kk) = sum(data(ik) .* wt(ik)) ./ sum(wt(ik));
		else
			centers(kk) = NaN;
		end
	end

	% assign any unused centers to the data point(s) furthest from centroids
	if redundant
		centers(isnan(centers)) = 0;
	elseif any(isnan(centers))
		printm('fixing unused centers %d', sum(isnan(centers)))
		while any(isnan(centers))
			cgood = col(centers(~isnan(centers)));
			index = quant1_index(data, cgood);
			dhat = cgood(index);
			iworst = imax(abs(data - dhat));
			knan = find(isnan(centers));
			centers(knan(1)) = data(iworst);
		end
	end
%	disp([iter centers])
	change = max(abs(centers - old));
	iter = iter + 1;
end
if iter == max_iter + 1
	warning 'max %d iterations reached'
end
if chat
	printf('%s: %d iterations', mfilename, iter)
end


% quant1_index()
% find index of nearest centroid.
% this version works even for complex data / centroids
%
function index = quant1_index(x, centers)
[dummy index] = min(abs(outer_sum(x, -centers)), [], 2);
%breaks = (centers(2:end) + centers(1:end-1)) / 2;
%index0 = 1 + sum(outer_sum(x, -breaks) > 0, 2);
%minmax(index-index0)


% quant1_rms()
% rms error between data and its quantized version (complex ok)
%
function rms = quant1_rms(x, centers)
index = quant1_index(x, centers);
rms = sqrt(mean(abs(x(:) - col(centers(index))).^2));


% lloyd_max_hist_test
% self test: compare this approach to matlab's lloyds routine
function lloyd_max_hist_test
rng(0)
x = [10*randn(10^4,1); 50 + 15*randn(10^4,1)];
%x = [10*rand(10^4,1); 50 + 15*rand(10^4,1)];

L = 5;
%c0 = 10*linspace(-1,1,L);
pn = jf_protected_names;
c0 = pn.prctile(x, 100*([1:L]-0.5)/L);

[nx cx] = hist(x, 100);

ch = highrate_centers(x, L);

if exist('lloyds') == 2
	tic
	[p1 c1] = lloyds(x, c0);
	t1 = toc;
else
	p1 = nan(L,1);
	c1 = nan(L,1);
	t1 = nan;
end

tic
c2 = lloyd_max_hist(x, c0, 50);
t2 = toc;

o = ones(L,1);
if im
	plot(cx, nx, '-', c0, 80*o, 'yo', ch, 60*o, 'yx', ...
		c1, 40*o, 'ys', c2, 20*o, 'y^')
	legend('hist', 'c0', 'highrate', 'lloyds', 'hist')
end

tic
c3 = lloyd_max_hist(x, c0);
t3 = toc;

if exist('quantiz') == 2
	[dum1 dum2 distor] = quantiz(x, p1, c1);
else
	distor = nan;
end

rms0 = quant1_rms(x, c0);
rms1 = quant1_rms(x, c1);
rms2 = quant1_rms(x, c2);
rms3 = quant1_rms(x, c3);
rmsh = quant1_rms(x, ch);
printf('lloyds matlab time=%g rms=%g distor=%g', t1, rms1, sqrt(distor))
printf('lloyd_max_hist time=%g rms=%g', t2, rms2)
printf('lloyd_max (no hist) time=%g rms=%g', t3, rms3)
printf('rms0=%g rmsh=%g', rms0, rmsh)
if im
	disp([c0, c1, c2, c3])
end
