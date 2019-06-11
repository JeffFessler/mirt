 function ob = Gdft(varargin)
%function ob = Gdft([args])
%| Construct Gdft object that computes samples of the DFT of a signal
%| with dimensions [(N)].  This is useful for "under-sampled" MRI.
%|
%| options (at least one of first two must be provided):
%|	'mask'	logical [(N)]	image-domain mask, usually: true(nx,ny)
%|	'samp'	logical [(N)]	which frequency-domain samples to return
%|				(default: all samples)
%|	'ifftshift'	0|1	apply ifftshift to image before DFT?
%|	'fftshift'	0|1	apply fftshift to spectrum after DFT?
%|				(both default to false)
%|
%| out
%|	ob	[nd np]		fatrix2 object
%|				nd = sum(samp(:)), np = sum(mask(:))
%|
%| See Gdft_test.m for example usage.
%|
%| Basically, you create a system matrix object by calling:
%|	A = Gdft( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically evaluate the DFT samples.
%| This is useful for iterative image reconstruction in MRI.
%|
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end

arg.mask = [];
arg.samp = [];
arg.ifftshift = false; % apply ifftshift to image before DFT?
arg.fftshift = false; % apply fftshift to spectrum after DFT?
arg.class = 'fatrix2';
%arg.class = 'Fatrix';
arg = vararg_pair(arg, varargin);

if isempty(arg.mask) && isempty(arg.samp), error 'must give mask or samp', end

if isempty(arg.mask)
	arg.mask = true(size(arg.samp));
elseif isempty(arg.samp)
	arg.samp = true(size(arg.mask));
end
if ~isequal(size(arg.mask), size(arg.samp)), fail 'size mismatch', end

if ~islogical(arg.mask), error 'mask must be logical', end
if ~islogical(arg.samp), error 'samp must be logical', end

switch arg.class
case 'Fatrix'

	arg.np = sum(arg.mask(:));
	arg.nd = sum(arg.samp(:));
	arg.ndim = ndims(arg.mask);
	if size(arg.mask, arg.ndim) == 1 % 1D case
		arg.ndim = arg.ndim - 1;
	end
	arg.dim = [arg.nd arg.np];
	ob = Fatrix(arg.dim, arg, ...
		'forw', @Gdft_forw_Fatrix, 'back', @Gdft_back_Fatrix, ...
		'gram', @Gdft_gram_Fatrix);

case 'fatrix2'
	forw = @(arg, x) arg.samp ...
		.* fftn_shifted(x, arg.ifftshift, arg.fftshift);
	back = @(arg, y) prod(size(arg.mask)) * arg.mask ...
		.* ifftn_shifted(y, arg.ifftshift, arg.fftshift);
	ob = fatrix2('arg', arg, 'imask', arg.mask, 'omask', arg.samp, ...
		'meth', {'plot', @Gdft_plot, ''}, ...
		'forw', forw, 'back', back);
otherwise
	fail('bad class')
end


% Gdft_abs()
% function out = Gdft_abs(ob)
% could be done, but probably useless


% Gdft_plot()
function out = Gdft_plot(arg)
im plc 1 2
im(1, arg.mask), title('mask')
im(2, arg.samp), title('samp')
out = [];


% Gdft_forw_Fatrix(): y = A * x
% in
%	x	[np L] or [(N) L]
% out
%	y	[M L]
%
function y = Gdft_forw_Fatrix(arg, x)

if size(x,1) == arg.np		% [np (L)]
	x = embed(x, arg.mask);	% [(N) (L)]
end

if isequal(size(x), size(arg.mask)) % [(N)]
	y = fftn_shifted(x, arg.ifftshift, arg.fftshift); % [(N)]
	y = y(arg.samp); % [M 1]
else
	xdim = size(x);
	NN = prod(size(arg.mask));
	LL = xdim(end);

	if arg.ndim == 1
		y = fft1_shifted(x, arg.ifftshift, arg.fftshift); % [(N) (L)]
		y = y(arg.samp,:); % [M (L)]

	else
		y = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(x, ll);
			y(:,ll) = col(fftn_shifted(tmp, ...
				arg.ifftshift, arg.fftshift));
		end
		y = y(arg.samp,:); % [M L]
	end
end


% Gdft_back_Fatrix(): x = A' * y
% in
%	y	[M L]
% out
%	x	[np L]
%
function x = Gdft_back_Fatrix(arg, y)

y = embed(y, arg.samp); % [(N) L]
ydim = size(y);
NN = prod(size(arg.mask));

if isequal(size(y), size(arg.samp)) % [(N)]
	x = ifftn_shifted(y, arg.ifftshift, arg.fftshift); % [(N)] adjoint
	x = col(x); % [*N]
else
	if arg.ndim == 1
		x = ifft1_shifted(y, arg.ifftshift, arg.fftshift); % [(N) L]
	else
		LL = ydim(end);
		x = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(y,ll);
			x(:,ll) = col(ifftn_shifted(tmp, ...
				arg.ifftshift, arg.fftshift));
		end
	end
end

% note the "NN" factor that is needed to make it an adjoint operation:
x = NN * x(arg.mask(:),:); % [np *L]

%x = reshape(x, [Ns Ld]); % [np (L)] % not needed if L can be scalar only!


% fft1_shifted()
function y = fft1_shifted(x, is_ifftshift, is_fftshift)
if is_ifftshift
	x = ifftshift(x, 1);
end
y = fft(x, [], 1);
if is_fftshift
	y = fftshift(y, 1);
end


% ifft1_shifted()
function x = ifft1_shifted(y, is_ifftshift, is_fftshift)
if is_fftshift
	y = fftshift(y, 1);
end
x = ifft(y, [], 1);
if is_ifftshift
	x = ifftshift(x, 1);
end


% fftn_shifted()
function y = fftn_shifted(x, is_ifftshift, is_fftshift)
if is_ifftshift
	x = ifftshift(x);
end
y = fftn(x);
if is_fftshift
	y = fftshift(y);
end


% ifftn_shifted()
function x = ifftn_shifted(y, is_ifftshift, is_fftshift)
if is_fftshift
	y = ifftshift(y);
end
x = ifftn(y);
if is_ifftshift
	x = fftshift(x);
end


% Gdft_gram_Fatrix()
function [T, reuse] = Gdft_gram_Fatrix(A, W, reuse)
if isempty(W)
	T = A' * A;
else
	T = A' * W * A;
end


function Gdft_test
%Nv = [10 6];
Nv = [7 6]; % test both odd and even sizes
%Nv = [4 5];
Gdft_test1(Nv)


% Gdft_test
function Gdft_test1(Nv)

rng(1)
samp = rand(Nv) > 0.3;
mask = true(Nv); mask(1:4) = false; % stress

cl_list = {'Fatrix', 'fatrix2'};
for ic = 1:numel(cl_list)
	cl = cl_list{ic};

	for ii = [1 0]
	for jj = [1 0]
		A = Gdft('mask', mask, 'samp', samp, 'class', cl, ...
			'ifftshift', ii, 'fftshift', jj);
		switch cl
		case 'Fatrix'
			Fatrix_test_basic(A, mask, 'complex', 1, 'halt', 0)
		case 'fatrix2'
			fatrix2_tests(A, 'complex', 1, 'halt', 0)
		end
		tolre = 2e-16; % Fatrix double
		tolre = 1e-7; % fatrix2 single
		test_adjoint(A, 'complex', 1, 'tolre', tolre);
	end
	end

	if 1 % test vs fft
		x = rand(Nv);
		x = x .* mask;
		y = fftn(x);
		y = y(samp(:));
		jf_equal(y, A * x(mask))

		x = ifftn(embed(y, samp)) * prod(Nv);
		x = x .* mask;
		jf_equal(x, embed(A' * y, mask))
	end

	if 1 % test gram
		wi = col(1:sum(samp(:)));

		switch cl
		case 'Fatrix'
			W = Gdiag(wi, 'class', 'Fatrix');
		case 'fatrix2'
			W = Gdiag(wi, 'mask', samp);
		end
		T = build_gram(A, W);

		% vector mode
		y1 = T * x(mask);
		y2 = (A' * (wi .* (A * x(mask))));
		jf_equal(y1, y2)

		% array mode
		y1 = T * x;
		switch cl
		case 'Fatrix'
			y2 = (A' * (wi .* (A * x)));
		case 'fatrix2'
			tmp = A * x;
			tmp = W * tmp;
			y2 = A' * tmp;
		end
		jf_equal(y1, y2)
	end
end
