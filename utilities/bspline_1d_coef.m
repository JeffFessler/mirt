 function coef = bspline_1d_coef(fn, varargin)
%function coef = bspline_1d_coef(fn, varargin)
%|
%| Given f[n], uniformly spaced samples of a 1D signal f(t),
%| compute the bspline coefficients so that the following model interpolates:
%| f(t) = \sum_{k=0}^{N-1} coef[k] b(t - k)
%|
%| in
%|	fn	[N (L)]		1d signal(s) of length N (columns)
%|
%| option
%|	order	default: 3
%|	ending	end / boundary conditions: mirror / periodic (default) / zero
%|
%| out
%|	coef	[N (L)]	bspline coefficients, for use in bspline_1d_interp()
%|
%| Caution: the 'zero' end conditions do not quite interpolate the edge samples.
%| That would require returning more coefficients than signal samples, which is
%| not worth the trouble since 'mirror' and 'periodic' are more common.
%|
%| Copyright 2005-12-7, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(fn, 'test'), bspline_1d_coef_test, return, end

arg.order = 3;
arg.ending = 'periodic';
arg.n0 = 10; % boundary conditions, more than enough for single precision.
arg.mex = 1; % default is to try mex version

arg = vararg_pair(arg, varargin);

dims = size(fn);

if dims(1) == 1 % handle row vector case
	if ndims(fn) > 2, fail 'unsupported', end
	N = dims(2);
else
	N = dims(1);
end

fn = reshapee(fn, N, []); % [N *L]

if arg.mex
	ii = strvcat('mirror', 'periodic', 'zero');
	ii = strmatch(arg.ending, ii, 'exact');
	if numel(ii) ~= 1, error 'unknown end conditions', end
	try
		coef = jf_mex('bspline,coef,1d', single(fn), ...
			int32(arg.order), int32(ii));
		coef = reshape(coef, dims); % [N (L)]
		return
	catch
		persistent warned
		if isempty(warned)
			warned = 1;
			warn 'mex version failed; reverting to matlab version'
		end
	end
end

switch arg.order
case 1
	coef = fn;
case 2
	coef = bspline2_1d_coef(fn, arg.n0, arg.ending);
case 3
	coef = bspline3_1d_coef(fn, arg.n0, arg.ending);
otherwise
	fail('order %d not done', arg.order)
end

coef = reshape(coef, dims); % [N (L)]


% bspline2_1d_coef()
% 1d quadratic bspline, various end conditions
% fn and coef are [N,L]
function coef = bspline2_1d_coef(fn, n0, ending)

p = -3 + sqrt(8); % pole location
fn = 8 * fn; % 1/a = 1/b(1)
coef = bspline_23_1d_coef(fn, n0, ending, p);


% bspline3_1d_coef()
% 1d cubic bspline, various end conditions
% fn and coef are [N,L]
function coef = bspline3_1d_coef(fn, n0, ending)

p = -2 + sqrt(3); % pole location
fn = 6 * fn; % 1/a = 1/b(1)
coef = bspline_23_1d_coef(fn, n0, ending, p);


% bspline_23_1d_coef()
% for 1d quadratic or cubic bspline, various end conditions
% fn and coef are [N,L]
function coef = bspline_23_1d_coef(fn, n0, ending, p)

N = size(fn,1);
coef = zeros(size(fn));

n0 = min(n0, N);
ps = p .^ [1:n0];

% initial condition for forward
switch ending
case 'mirror'
	coef(1,:) = ps/p * fn(1:n0,:);
case 'periodic'
	coef(1,:) = fn(1,:) + fliplr(ps) * fn([(N-n0+1):N],:);
case 'zero'
	coef(1,:) = fn(1,:);
otherwise
	fail('unknown end condition "%s"', ending)
end

% causal iir filter
for n=2:N
	coef(n,:) = fn(n,:) + p * coef(n-1,:);
end

% initial condition for anti-causal
if streq(ending, 'mirror')
%	coef(N,:) = -fliplr(ps) * coef([(N-n0+1):N],:); % my way, not good
	coef(N,:) = -p/(1-p^2) * ( coef(N,:) + p * coef(N-1,:) ); % unser way
elseif streq(ending, 'periodic')
	coef(N,:) = -p * coef(N,:) - p * ps * coef(1:n0,:);
else % zero
	coef(N,:) = -p * coef(N,:);
end

% anti-causal iir filter
for n=N-1:-1:1
	coef(n,:) = p * ( coef(n+1,:) - coef(n,:) );
end


% fft approach, for validating periodic and mirror cases
function coef = bspline2_1d_coef_fft(fn, ending)
N = size(fn,1);
M = size(fn,2);

switch ending
case 'mirror'
	hn = zeros(2*N-2,1);
	hn([1 2 2*N-2]) = [6 1 1]'/8;
	coef = ifft(fft([fn; flipud(fn(2:N-1,:))]) ./ ...
		repmat(fft(hn), [1 M]));
	coef = reale(coef(1:N,:));
case 'periodic'
	hn = zeros(N,1);
	hn([1 2 N]) = [6 1 1]'/8;
	coef = ifft(fft(fn) ./ repmat(fft(hn), [1 M]));
	coef = reale(coef);
case 'zero'
	coef = zeros(size(fn)); % fake
otherwise
	fail 'bug'
end


% fft approach, for validating periodic and mirror cases
function coef = bspline3_1d_coef_fft(fn, ending)
N = size(fn,1);
M = size(fn,2);

switch ending
case 'mirror'
	hn = zeros(2*N-2,1);
	hn([1 2 2*N-2]) = [4 1 1]'/6;
	coef = ifft(fft([fn; flipud(fn(2:N-1,:))]) ./ ...
		repmat(fft(hn), [1 M]));
	coef = reale(coef(1:N,:));
case 'periodic'
	hn = zeros(N,1);
	hn([1 2 N]) = [4 1 1]'/6;
	coef = ifft(fft(fn) ./ repmat(fft(hn), [1 M]));
	coef = reale(coef);
case 'zero'
	coef = zeros(size(fn)); % fake
otherwise
	fail 'bug'
end


% bspline3_1d_coef_exact()
% exact matrix approach based on cbanal.m from kybic, for validating
function coef = bspline3_1d_coef_exact(y, ending)
N = size(y,1);
A = zeros(N,N);

if N == 1
	coef = 1.5 * y; % 6/4
	return
end

if streq(ending, 'mirror')
	coef = bspline3_1d_coef_exact([y; flipud(y(2:N-1,:))], 'periodic');
	coef = coef(1:N,:);
	return
end

A(1,1:2) = [4 1]/6;
A(N,N-1:N) = [1 4]/6;
if streq(ending, 'mirror')
	A(1,2) = 2/6;
	A(N,N-1) = 2/6;
elseif streq(ending, 'periodic')
	A(1,N) = 1/6;
	A(N,1) = 1/6;
end
for i=2:N-1,
	A(i,i-1:i+1) = [1 4 1]/6;
end
coef = A \ y;


% bspline_1d_coef_test()
% test periodic case by comparing to fft and to matrix "exact" method
function bspline_1d_coef_test

N = 2^5;
%fn = zeros(N,1); fn(3) = 1;
fn = eye(N);

orders = [1 2 3];
endings = {'mirror', 'periodic', 'zero'};
for io=1:numel(orders)
	order = orders(io);

	for ie=1:numel(endings)
		ending = endings{ie};
		printm([ending ' ' num2str(order)])

		arg = {'order', order, 'ending', ending};
		if 1 % check mex version
%			cpu etic
			coef = bspline_1d_coef(fn, arg{:}, 'mex', 0);
%			cpu etoc 'matlab'
%			cpu etic
			cmex = bspline_1d_coef(fn, arg{:}, 'mex', 1);
%			cpu etoc 'mex'
			equivs(coef, cmex) % max_percent_diff(coef, cmex)
			if im
				im plc 1 3
				im(1, coef)
				im(2, cmex)
				im(3, cmex-coef)
				xlabelf('order %d, ending %s', order, ending)
				drawnow
			end
		end

		if has_mex_jf % check adjoint
			A = bspline_1d_coef(eye(N), arg{:}, 'mex', 1);
			B = jf_mex('bspline,coef,adj,1d', single(eye(N)), ...
				int32(order), int32(ie));
			if max(col(abs(A-B'))) > 1e-5, error 'adjoint', end
		end

		if order ~= 1 && im, prompt, end
		if order ~= 3, continue, end

		cfft = bspline3_1d_coef_fft(fn, ending);
		exact = bspline3_1d_coef_exact(fn, ending);

		if ~streq(ending, 'zero')
			max_percent_diff(exact, cfft) % tiny, except for 'zero'
		end
		max_percent_diff(exact, coef) % todo: 7% error for 'zero' !?

		if im
			im plc 2 3
			im(4, fn, 'f[n]'), cbar
			xlabelf('order %d', order)
			im(1, exact, 'exact'), cbar
			xlabelf('ending %s', ending)
			im(2, coef, 'coef'), cbar
			im(5, coef-exact, 'error'), cbar
			if ~streq(ending, 'zero')
				im(3, cfft, 'coef fft'), cbar
				im(6, cfft-exact, 'error'), cbar
			end
		prompt
		end
	end
end
