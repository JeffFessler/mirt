 function ft = bspline_1d_interp(fn, ti, varargin)
%function ft = bspline_1d_interp(fn, ti, varargin)
%|
%| given f[n], the "unit spaced" samples of a continuous signal f(t),
%| perform b-spline interpolation at sample locations t_i using the model:
%|	f(t) = \sum_{k=0}^{N-1} coef_k b(t - k)
%| where the coefficients are computed internally using bspline_1d_coef.m
%|
%| in
%|	fn	[N (L)]			1d signal(s) samples
%|	ti	[M 1] or [M (L)]	desired sample points (unitless)
%|
%| option
%|	order		default: 3
%|	ending		end / boundary conditions: mirror / periodic / zero
%|				default: 'periodic'
%|	mex		0 to disable mex (default: 1)
%|	ob		1 to create Fatrix or fatrix2 object (default: 0)
%|	'class'		'fatrix2' (default) or 'Fatrix' (obsolete)
%|
%| out
%|	ft	[M (L)]			f(t_i) interpolated signal values
%|
%| Copyright 2005-12-7, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fn, 'test'), bspline_1d_interp_test, return, end
if nargin < 2, ir_usage, end

arg.order = 3;
arg.ending = 'periodic';
arg.mex = 1; % default is to try mex
arg.ob = false; % 1 to create fatrix2 or Fatrix
arg.class = 'fatrix2';

arg = vararg_pair(arg, varargin);

arg.ti = ti;
arg.N = size(fn,1);

ii = strvcat('mirror', 'periodic', 'zero');
arg.ie = strmatch(arg.ending, ii, 'exact');
if numel(arg.ie) ~= 1, error 'unknown end conditions', end

if arg.ob
	ft = bspline_1d_interp_ob(arg);
return
end

ft = bspline_1d_interp_arg(arg, fn);


% bspline_1d_interp_ob()
function ob = bspline_1d_interp_ob(arg)
if size(arg.ti,2) == 1 % usual 1D case
	odim = numel(arg.ti);
	idim = arg.N;
	does_many = true;
else % trickier case
	odim = size(arg.ti);
	idim = [arg.N size(arg.ti,2)];
	does_many = false; % unsure, so play it safe
end

switch arg.class
case 'Fatrix'
	dim = [prod(odim) prod(idim)];
	ob = Fatrix(dim, arg, 'caller', 'bspline_1d_interp', ...
		'forw', @bspline_1d_interp_arg, 'back', @bspline_1d_interp_adj);

case 'fatrix2'
	ob = fatrix2('arg', arg, 'does_many', does_many, ...
		'odim', odim, 'idim', idim, ...
		'forw', @bspline_1d_interp_arg, 'back', @bspline_1d_interp_adj);

otherwise
	fail 'class'
end


% bspline_1d_interp_adj()
function fn = bspline_1d_interp_adj(arg, ft)

if arg.mex
	try
		fn = jf_mex('bspline,interp,adj,1d', ...
			single(ft), single(arg.ti), int32(arg.N), ...
			int32(arg.order), int32(arg.ie));
		return
	catch
		persistent warned
		if isempty(warned)
			warned = 1;
			warn 'mex failed; revert to non-mex'
		end
	end
end

fail 'non-mex adjoint not done' % if added later, then add test too!


% bspline_1d_interp_arg()
function ft = bspline_1d_interp_arg(arg, fn)

if arg.mex
	try
		ft = jf_mex('bspline,interp,1d', single(fn), single(arg.ti), ...
			int32(arg.order), int32(arg.ie));
		return
	catch
		persistent warned
		if isempty(warned)
			warned = 1;
			warn 'mex failed; revert to non-mex'
		end
	end
end

ck = bspline_1d_coef(fn, 'order', arg.order, 'ending', arg.ending);
ft = bspline_1d_synth(ck, arg.ti, 'order', arg.order, 'ending', arg.ending);


% bspline_1d_interp_test
% test
function bspline_1d_interp_test

N = 12;
ff = @(t) 1 + cos(2*pi*t/6);
n = [0:N-1]';
fn = ff(n);

ti = linspace(-4, N+4, 41*(N+8)+1)';
ft = ff(ti);
if 0 % test multi
	ti = [ti flipud(ti)];
	fn = [fn 2+flipud(fn)];
end

%fn = single(fn);
%ti = single(ti);

orders = [1 2 3];
endings = {'mirror', 'periodic', 'zero'};

fx = cell(numel(orders), numel(endings));
fm = cell(numel(orders), numel(endings));
for io = 1:numel(orders)
	order = orders(io);
	for ie = 1:numel(endings)
		ending = endings{ie};
		arg = {'order', order, 'ending', ending};
		fx{io,ie} = bspline_1d_interp(fn, ti, arg{:}, 'mex', 1);
		fm{io,ie} = bspline_1d_interp(fn, ti, arg{:}, 'mex', 0);
%		max_percent_diff(fx{io,ie}, fm{io,ie}) % mex vs mat
		equivs(fx{io,ie}, fm{io,ie}) % mex vs mat
		if 0 && io==1 && ie==1
			plot(n, fn, 'ro', ti, ft, 'y-', ...
				ti, fx{io,ie}, 'g', ti, fm{io,ie}, 'c')
			keyboard
		end

		A = bspline_1d_interp(single(eye(N)), ti, arg{:});

		% adjoint done only for mex version, so tested only that way
		if has_mex_jf % test mex adjoint
			B = jf_mex('bspline,interp,adj,1d', ...
				single(eye(numel(ti))), ...
				single(ti), int32(N), int32(order), int32(ie));
			if max(col(abs(A-B'))) > 1e-5, error 'adjoint', end
		end

		% test 'ob' version
		classes = {'Fatrix', 'fatrix2'};
		for ic = 1:numel(classes)
			B = bspline_1d_interp(nan(N,1), ti, arg{:}, ...
				'ob', 1, 'class', classes{ic});
			B = B * eye(N);
			equivs(A, B)
		end

		if has_mex_jf % test adjoint of 'ob' version
			A = bspline_1d_interp(ones(N,1), ti, arg{:}, 'ob', 1);
			B = A' * eye(numel(ti));
			A = A * eye(N);
			if max(col(abs(A-B'))) > 1e-5, error 'adjoint', end
		end
	end

	if im
		subplot(310+io)
		plot(n, fn, 'ro', ti, ft, 'b-', ...
			ti, fm{io,1}, 'm-.', ...
			ti, fm{io,2}, 'y--', ...
			ti, fm{io,3}, 'g--')
		legend('sample', 'f(t)', 'm', 'p', 'z', 'location', 'east')
		ylabelf('order = %d', io), ytick([0 2]), axis([-4 21 -0.1 2.3])
	end
end

%max_percent_diff(ft, fm{2,2})
