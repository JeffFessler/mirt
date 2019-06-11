 function ob = Gdsft(om, Nd, varargin)
%function ob = Gdsft(om, Nd, varargin)
%| Construct Gdsft object, which computes (nonuniform) FT samples
%| of signals with dimensions [(Nd)] exactly, e.g., for testing Gnufft
%| For faster computation, use Gnufft instead.
%| See Gdsft_test.m for example usage.
%|
%| in
%|	om	[M D]		frequency locations (radians / sample)
%|	Nd	[1 D]		signal dimensions
%|
%| options
%|	mask	[(Nd)]		logical support array
%|	n_shift	[1 D]		see nufft_init
%|	nthread	[]		# of processor threads
%|	class	''		fatrix2 (default) or Fatrix or 'exact'
%|	use_mex	0|1		use jf_mex ? (default: has_mex_jf)
%|
%| out
%|	ob	[M np]		object
%|
%| Copyright 2005-7-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(om, 'test'), Gdsft_test, return, end
if nargin < 2, ir_usage, end

% defaults
arg.om_t = double(om)'; % [D M], as required by dtft_mex
arg.Nd = Nd;
arg.mask = [];
arg.n_shift = '';
arg.nthread = 1;
arg.class = 'fatrix2';
arg.use_mex = has_mex_jf;
arg.show = false; % used in Gdsft_gram()

% options
arg = vararg_pair(arg, varargin);
arg.ndim = numel(arg.Nd);

if isempty(arg.n_shift)
	arg.n_shift = zeros(1, arg.ndim);
end

if size(arg.om_t,1) ~= arg.ndim
	error 'dimension mismatch'
end

switch arg.class
case 'Fatrix'
	ob = Gdsft_Fatrix(arg);

case 'fatrix2'
	ob = Gdsft_fatrix2(arg);

case 'exact'
	ob = Gdsft_fatrix2(arg);
	ob = full(ob); % trick

otherwise
	fail 'class'
end


% Gdsft_fatrix2()
function ob = Gdsft_fatrix2(arg)
if ~arg.use_mex
	forw = @(arg, x) dtft(x, arg.om_t', 'n_shift', arg.n_shift);
	back = @(arg, y) fatrix2_maskit(arg.mask, ...
		dtft_adj(y, arg.om_t', arg.Nd, arg.n_shift));

elseif any(arg.n_shift ~= 0)
	arg.phasor = exp(1i * (arg.om_t' * arg.n_shift(:))); % [M 1]
	forw = @(arg, x) arg.phasor .* ...
		cast(jf_mex('dtft,forward', arg.om_t, ...
			double(x), int32(arg.nthread)), class(x));
	back = @(arg, y) fatrix2_maskit(arg.mask, ...
		cast(jf_mex('dtft,adjoint', arg.om_t, ...
			complexify(double(conj(arg.phasor) .* y)), ...
			int32(arg.Nd), int32(arg.nthread)), class(y)));
else
	forw = @(arg, x) ...
		cast(jf_mex('dtft,forward', arg.om_t, ...
			double(x), int32(arg.nthread)), class(x));
	back = @(arg, y) fatrix2_maskit(arg.mask, ...
		cast(jf_mex('dtft,adjoint', arg.om_t, ...
			complexify(double(y)), ...
			int32(arg.Nd), int32(arg.nthread)), class(y)));
end
odim = size(arg.om_t,2);
ob = fatrix2('arg', arg, 'imask', arg.mask, ...
	'idim', arg.Nd, 'odim', odim, ...
	'gram', @Gdsft_gram, 'forw', forw, 'back', back);


% Gdsft_Fatrix()
function ob = Gdsft_Fatrix(arg)

arg = Gdsft_init_Fatrix(arg); % initialize

if isempty(arg.mask)
	arg.mask = true([arg.Nd 1]); % [(Nd)]
end
arg.np = sum(arg.mask(:));
arg.dim = [size(arg.om_t,2) arg.np]; % [M np]
ob = Fatrix(arg.dim, arg, ...
	'gram', @Gdsft_gram_fail,  ...
	'forw', @Gdsft_forw_Fatrix, 'back', @Gdsft_back_Fatrix);


% Gdsft_init_Fatrix()
function arg = Gdsft_init_Fatrix(arg)

if ~isempty(arg.n_shift)
% fix: - ?
	arg.phasor = exp(1i * (arg.om_t' * arg.n_shift(:))); % [M 1]
	arg.phasor = diag_sp(arg.phasor); % trick: to handle multiples
else
	arg.phasor = 1;
end


% Gdsft_forw_Fatrix(): y = A * x
% in
%	x	[np L] or [(Nd) L]
% out
%	y	[M L]
%
function y = Gdsft_forw_Fatrix(arg, x)

if size(x,1) == arg.np
	x = embed(x, arg.mask);	% [(Nd) (L)]
end

if arg.use_mex
	y = jf_mex('dtft,forward', arg.om_t, double(x), int32(arg.nthread));
else
	y = dtft(x, arg.om_t', 'n_shift', arg.n_shift);
end
y = arg.phasor * y;


% Gdsft_back_Fatrix(): x = A' * y
% in
%	y	[M L]
% out
%	x	[np L]
%
function x = Gdsft_back_Fatrix(arg, y)

y = arg.phasor' * y;
if isreal(y)
	y = complexify(y);
end
if arg.use_mex
	x = jf_mex('dtft,adjoint', arg.om_t, y, int32(arg.Nd), int32(arg.nthread));
else
	x = dtft_adj(y, arg.om_t', arg.Nd, arg.n_shift);
end
x = reshape(x, prod(arg.Nd), []);
x = x(arg.mask(:),:);


% Gdsft_gram()
function [T, reuse] = Gdsft_gram_fail(ob, W, reuse)
% T = dsft_gram(ob, W);
error 'not done'
