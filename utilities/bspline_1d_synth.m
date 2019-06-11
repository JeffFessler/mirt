 function ft = bspline_1d_synth(coef, ti, varargin)
%function ft = bspline_1d_synth(coef, ti, varargin)
%|
%| Given c[n], the b-spline coefficients for the following model:
%| f(t) = \sum_{k=0}^{N-1} coef_k b(t - k),
%| and given sorted sample locations ti, synthesize signal samples f(t_i).
%|
%| in
%|	coef	[N (L)]			bspline coefficients for k=0,...,N-1
%|					usually from bspline_1d_coef
%|	ti	[M 1] or [M (L)]	desired sample points (unitless)
%|
%| option
%|	order	1|2|3	b-spline order, default: 3
%|	ending		end / boundary conditions: mirror / periodic / zero
%|			default: 'periodic'
%|	mex	0|1	use mex version? default: 1
%|
%| out
%|	ft	[M (L)]			f(t_i) interpolated signal values
%|
%| Copyright 2005-12-7, Jeff Fessler, University of Michigan

if nargin == 1 && streq(coef, 'test'), im clf, bspline_1d_synth_test, return, end
if nargin == 1 && streq(coef, 'plot'), im clf, bspline_1d_synth_plot, return, end
if nargin < 2, ir_usage, end

arg.order = 3;
arg.ending = 'periodic';
arg.mex = 1; % default is to try mex version

arg = vararg_pair(arg, varargin);

if arg.mex
	try
		ft = jf_mex('bspline,synth,1d', single(coef), single(ti), ...
			int32(arg.order), ...
			bspline_1d_synth_mex_ending(arg.ending));
		return
	catch
		persistent warned
		if isempty(warned)
			warned = 1;
			warn 'mex version failed; reverting to matlab version'
		end
	end
end

dims = size(coef);
if dims(1) == 1 % handle row vector case
	if ndims(coef) > 2, fail 'unsupported', end
	N = dims(2);
	if sum(size(ti) > 1) > 1, fail 'unsupported', end
	ti = ti(:); % [M 1]
	dims = [1 numel(ti)]; % output dim is also a row vector
else
	N = dims(1);
	dims(1) = size(ti,1); % output dim is [M (L)]
end

coef = reshape(coef, N, []); % [N *L]
nc = ncol(coef);

[coef ti] = bspline_1d_setup(coef, ti, arg.order, arg.ending);

if (ncol(ti) ~= 1 && ncol(ti) ~= ncol(coef)), error 'ti size mismatch', end

ft = zeros(size(ti,1), size(coef,2));
switch arg.order
case 1
	for ic = 1:nc
		it = 1 + (ic <= ncol(ti)) * (ic-1); % 1 or ic
		ft(:,ic) = bspline1_1d_synth(coef(:,ic), ti(:,it));
	end

case 2
	for ic = 1:nc
		it = 1 + (ic <= ncol(ti)) * (ic-1); % 1 or ic
		ft(:,ic) = bspline2_1d_synth(coef(:,ic), ti(:,it));
	end

case 3
	for ic = 1:ncol(coef)
		it = 1 + (ic <= ncol(ti)) * (ic-1); % 1 or ic
		ft(:,ic) = bspline3_1d_synth(coef(:,ic), ti(:,it));
	end

otherwise
	fail('order %d not done', arg.order)
end

ft = reshape(ft, dims); % [N (L)]


% bspline_1d_setup()
% trick: preprocess locations and pad coefficients to avoid conditionals!
% and to handle boundary conditions cleanly.
function [coef, ti] = bspline_1d_setup(coef, ti, order, ending)

[N L] = size(coef);

if streq(ending, 'mirror')
	% map ti into [0,N-1] using mirror conditions
	ti = mod(ti, 2*(N-1)); ti(ti > N-1) = 2*(N-1)-ti(ti > N-1);
	switch order
	case 1
		coef = [coef; zeros(1,L)];
	case {2, 3}
		coef = [coef(2,:); coef; coef(N-1,:); zeros(1,L)];
		ti = ti + 1;
	otherwise
		fail('order %d not done', order)
	end

elseif streq(ending, 'periodic')
	ti = mod(ti, N); % now ti is in [0,N)
	switch order
	case 1
		coef = [coef; coef(1,:)];
	case {2, 3}
		coef = [coef(N,:); coef; coef(1:2,:)];
		ti = ti + 1;
	otherwise
		fail('order %d not done', order)
	end

elseif streq(ending, 'zero')
	switch order
	case 1
		ti = max(ti, -1);
		ti = min(ti, N);
		coef = [zeros(1,L); coef; zeros(2,L)];
		ti = ti + 1; % [0,N+1]
	case 2
		ti = max(ti, -1.5);
		ti = min(ti, N+0.5);
		coef = [zeros(2,L); coef; zeros(3,L)];
		ti = ti + 2; % [0.5,N+2.5]
	case 3
		ti = max(ti, -2);
		ti = min(ti, N+1);
		coef = [zeros(3,L); coef; zeros(4,L)];
		ti = ti + 3; % [1,N+4]
	otherwise
		fail('order %d not done', order)
	end

end


% bspline1_1d_synth()
% 1d linear bspline
% coef [N L]	b-spline coefficients. trick: with padding at head and tail!
% ti [M 1]	trick: must be in range [0,N-2]
% ft [M L]
function ft = bspline1_1d_synth(coef, ti)

n0 = floor(ti);
ft = coef(1+n0) .* (1 - (ti - n0)) + coef(2+n0) .* (ti - n0);


% bspline2_1d_synth()
% 1d quadratic bspline
% coef [N L]	b-spline coefficients. trick: with padding at head and tail!
% ti [M 1]	trick: must be in range [1/2,N-7/2)
% ft [M L]
function ft = bspline2_1d_synth(coef, ti)

N = size(coef,1);
M = length(ti);
n0 = floor(ti - 1/2);

b2f0 = @(t) 3/4 - t.^2;
b2f1 = @(t) (abs(t) - 3/2).^2 / 2;

ft =	coef(1+n0) .* b2f1(ti - (n0)) + ...
	coef(2+n0) .* b2f0(ti - (n0+1)) + ...
	coef(3+n0) .* b2f1(ti - (n0+2));


% bspline3_1d_synth()
% 1d cubic bspline
% coef [N L], ti [M 1]
% ft [M L]
% coef [N L]	b-spline coefficients. trick: with padding at head and tail!
% ti [M 1]	trick: must be in range [1,N-4)
% ft [M L]
function ft = bspline3_1d_synth(coef, ti)

N = size(coef,1);
M = length(ti);
n0 = floor(ti-1);
%if any(n0 < 1) || any(3+n0 > N), minmax(n0), N, keyboard, error, end

ft =	coef(1+n0) .* b3f0(ti - (n0)) + ...
	coef(2+n0) .* b3f1(ti - (n0+1)) + ...
	coef(3+n0) .* b3f2(ti - (n0+2)) + ...
	coef(4+n0) .* b3f3(ti - (n0+3));


% bspline3_1d_synth_old()
% 1d cubic bspline.  this version does not require padding, but is slow
% coef [N L], ti [M 1]
% ft [M L]
function ft = bspline3_1d_synth_old(coef, ti, is_periodic)

N = size(coef,1);
M = length(ti);
ft = zeros(M, ncol(coef));
n0 = floor(ti);

i0 = (n0 >= 1) & (n0 < N+1);
i1 = (n0 >= 0) & (n0 < N);
i2 = (n0 >= -1) & (n0 < N-1);
i3 = (n0 >= -2) & (n0 < N-2);

for ic = 1:ncol(coef)
	ft(i0,ic) =		coef(0+n0(i0),ic) .* b3f0(ti(i0) - (n0(i0)-1));
	ft(i1,ic) = ft(i1,ic) + coef(1+n0(i1),ic) .* b3f1(ti(i1) - n0(i1));
	ft(i2,ic) = ft(i2,ic) + coef(2+n0(i2),ic) .* b3f2(ti(i2) - (n0(i2)+1));
	ft(i3,ic) = ft(i3,ic) + coef(3+n0(i3),ic) .* b3f3(ti(i3) - (n0(i3)+2));

	if is_periodic
		ie = (n0 == 0);
		ft(ie,ic) = ft(ie,ic) + coef(N,ic) .* b3f0(ti(ie) + 1);
		ie = (n0 == N-1);
		ft(ie,ic) = ft(ie,ic) + coef(1,ic) .* b3f2(ti(ie) - N);
		ft(ie,ic) = ft(ie,ic) + coef(2,ic) .* b3f3(ti(ie) - (N+1));
		ie = (n0 == N-2);
		ft(ie,ic) = ft(ie,ic) + coef(1,ic) .* b3f3(ti(ie) - N);
	end
end

% piece of bspline function for each offset
%b3f = @(t) 1/6*(2 - abs(t)).^3 .* (1 < abs(t) & abs(t) < 2) ...
%	+ (2/3 - t.^2 + 1/2 * abs(t).^3) .* (abs(t) <= 1);

function out = b3f0(t), out = (2 - t).^3 / 6;
%function out = b3f1(t), out = 2/3 - t.^2 + 1/2 * t.^3;
function out = b3f1(t), out = 2/3 - t.^2 .* (1 - t/2);
%function out = b3f2(t), out = 2/3 - t.^2 - 1/2 * t.^3;
function out = b3f2(t), out = 2/3 - t.^2 .* (1 + t/2);
function out = b3f3(t), out = (2 + t).^3 / 6;


% bspline_1d_synth_mex_ending()
% mex file expects an integer for the ending option
% 1=mirror 2=periodic 3=zero
function ending = bspline_1d_synth_mex_ending(ending)
list = strvcat('mirror', 'periodic', 'zero');
ii = strmatch(ending, list, 'exact');
if length(ii) ~= 1, error 'unknown end conditions', end
ending = int32(ii);


% test
function bspline_1d_synth_test

N = 12;
fn = eye(N);
ti = linspace(-3,N+3,41*(N+6)+1)';
M = length(ti);

orders = [1 2 3];
endings = {'mirror', 'periodic', 'zero'};
for io=1:length(orders)
	order = orders(io);

	for ie=1:length(endings)
		ending = endings{ie};
		printm([ending ' ' num2str(order)])

		arg = {'order', order, 'ending', ending};
%		coef = bspline_1d_coef(fn, arg{:});
		coef = eye(N);

		if 1
			fmat = bspline_1d_synth(coef, ti, arg{:}, 'mex', 0);
			fmex = bspline_1d_synth(coef, ti, arg{:}, 'mex', 1);
			if im
				clf, plot(ti, fmat(:,1), 'g-', ti, fmex(:,1), 'b-')
			end
%			prompt
			equivs(fmex, fmat) % max_percent_diff(fmex,fmat)
		end

		if has_mex_jf % test adjoint
			A = bspline_1d_synth(eye(N), ti, arg{:}, 'mex', 1);
			B = jf_mex('bspline,synth,adj,1d', single(eye(M)), ...
				single(ti), int32(N), int32(order), ...
				bspline_1d_synth_mex_ending(ending))';

			if max(col(abs(A-B))) > 1e-5, error 'adjoint', end
		end
	end
end

if im
	clf
	pl = @(i) subplot(310+i);
	n = 0:N-1;
	pl(1), plot(n, fn, '+', n, coef, 'o', ti, fmat, '-')
	axis([-3 N+3 0 1.1]), xtick([-2 0 N-1 N+1])
	pl(2), plot(n, sum(fn,2), '+', ti, sum(fmat,2), '-')
	axis([-3 N+3 0 1.1]), xtick([-2 0 N-1 N+1])
	pl(3), bspline_1d_synth_plot
end


% bspline_1d_synth_plot
function bspline_1d_synth_plot
ti = linspace(-1,5,601)';
for order = 1:3
	bn(:,order) = bspline_1d_synth([0 0 1 0 0]', ti, 'order', order, ...
		'mex', 0, 'ending', 'zero');
end
plot(ti-2, bn)
ytick([0 1/8 1/6 2/3 3/4 1])
set(gca, 'yticklabel', strvcat({'0', '1/8', '1/6', '2/3', '3/4', '1'}))
legend('1', '2', '3')
grid
