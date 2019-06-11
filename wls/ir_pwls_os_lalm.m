 function [xs, info] = ir_pwls_os_lalm(x, Ab, yi, R, varargin)
%function [xs, info] = ir_pwls_os_lalm(x, Ab, yi, R, [options])
%|
%| (This is now superceded by new method in ir_pwls_os_rlalm.m)
%|
%| penalized weighted least squares estimation / image reconstruction
%| using linearized augmented lagrangian method with
%| (optionally relaxed) ordered subsets. (OS-LALM)
%|
%| See Feb. 2015 IEEE T-MI by Hung Nien & J A Fessler
%| "Fast X-ray CT image reconstruction using
%|	the linearized augmented Lagrangian method with ordered subsets"
%|	http://dx.doi.org/10.1109/TMI.2014.2358499
%|
%| cost(x) = (y-Ax)' W (y-Ax) / 2 + R(x)
%|
%| in
%|	x	[np 1]		initial estimate
%|	Ab	[nd np]		Gblock object (needs abs(Ab) method)
%|				or sparse matrix (implies nsubset=1)
%|	yi	[nb na]		measurements (noisy sinogram data)
%|	R	penalty		object (see Reg1.m), can be []
%|
%| option
%|	niter			# of iterations (default: 1)
%|	nfista			# of fista iterations (default: 1)
%|	wi	[nb na]		weighting sinogram (default: [] for uniform)
%|	pixmax	[1] or [2]	max pixel value, or [min max] (default [0 inf])
%|	denom	[np 1]		precomputed denominator
%|	aai	[nb na]		precomputed row sums of |Ab|
%|	relax0	[1] or [2]	relax0 or (relax0, relax_rate) (default 1)
%|	userfun	@		user-defined function handle (see default below)
%|					taking arguments (x, userarg{:})
%|	userarg	{}		user arguments to userfun (default {})
%|	chat
%|
%| out
%|	xs	[np niter]	iterates
%|	info	[niter 1]	time
%|
%| 2014-09-12, Hung Nien, based on pwls_os_sqs
%| 2014-09-15 tweaks by Jeff Fessler
%| 2015-08-16, Jeff Fessler, added support for Ab with negatives

if nargin == 1 && streq(x, 'test'), ir_pwls_os_lalm_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

% defaults
arg.niter = 1;
arg.nfista = 1;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.pixmax = inf;
arg.chat = false;
arg.wi = [];
arg.aai = [];
arg.rho0 = 1; % initial value for rho
arg.relax0 = 1;
arg.denom = [];
arg.scale_nblock = true; % traditional scaling
arg.update_even_if_denom_0 = true;
arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

Ab = block_op(Ab, 'ensure'); % make it a block object (if not already)
nblock = block_op(Ab, 'n');
starts = subset_start(nblock);

cpu etic

wi = arg.wi;
if isempty(wi)
	wi = ones(size(yi));
end
if isempty(arg.aai) && isempty(arg.denom)
	arg.aai = reshape(sum(abs(Ab)'), size(yi)); % a_i = sum_j |a_ij|
end

% check input sinogram sizes for OS
if (ndims(yi) ~= 2) || (size(yi,2) == 1 && nblock > 1)
	fail 'bad yi size'
end
if (ndims(wi) ~= 2) || (size(wi,2) == 1 && nblock > 1)
	fail 'bad wi size'
end

relax0 = arg.relax0(1);
if length(arg.relax0) == 1
	relax_rate = 0;
elseif length(arg.relax0) == 2
	relax_rate = arg.relax0(2);
else
	error relax
end

if length(arg.pixmax) == 2
	pixmin = arg.pixmax(1);
	pixmax = arg.pixmax(2);
elseif length(arg.pixmax) == 1
	pixmin = 0;
	pixmax = arg.pixmax;
else
	error pixmax
end

% likelihood denom, if not provided
denom = arg.denom;
if isempty(denom)
	denom = abs(Ab)' * col(arg.aai .* wi); % needs abs(Ab)
	arg = rmfield(arg, 'aai'); % clear to save memory - not needed below
end
if ~arg.update_even_if_denom_0
	% todo: this may not work for LALM because "denom" appears in numerator!
	denom(denom == 0) = inf; % trick: prevents pixels where denom=0 being updated
end

if isempty(R)
	pgrad = 0; % unregularized default
	Rdenom = 0;
end

[nb na] = size(yi);

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave), 'single');
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = single(x);
end

%info = zeros(niter,?); % do not initialize since size may change

% initialization
rho = arg.rho0;

iblock = nblock;
ia = iblock:nblock:na;
li = Ab{iblock} * x;
li = reshape(li, nb, length(ia));
resid = wi(:,ia) .* (li - yi(:,ia));
if arg.scale_nblock
	scale = nblock; % traditional way
else
	scale = na / numel(ia); % alternative - untested
end
zeta = scale * Ab{iblock}' * resid(:);

g = zeta;

ii = 1;

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	relax = relax0 / (1 + relax_rate * (iter-1));

	% loop over subsets
	for iset = 1:nblock
		s = rho * zeta + (1-rho) * g;

		% inner denoising problem
		z = x; v = x;
		zold = z; told = 1;
		% fast iterative shrinkage/thresholding algorithm (FISTA)
		for ifista = 1:arg.nfista
			sigma = rho * denom .* (v-x) + s;
			if ~isempty(R)
				pgrad = R.cgrad(R, v);
				Rdenom = R.denom(R, v);
			end

			num = sigma + pgrad;
			den = rho * denom + Rdenom;

			z = v - relax * num ./ den;
			z = max(z, pixmin);
			z = min(z, pixmax);

			% adaptive restart for FISTA
			if (v-z)' * (z-zold) > 0
				t = 1;
				v = z;
			else
				t = (1 + sqrt(1 + 4 * told^2)) / 2;
				v = z + (told-1) / t * (z-zold);
			end

			zold = z; told = t;
		end
		x = z;

		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Ab{iblock} * x;
		li = reshape(li, nb, length(ia));
		resid = wi(:,ia) .* (li - yi(:,ia));

		if arg.scale_nblock
			scale = nblock; % traditional way
		else
			scale = na / numel(ia); % alternative - untested
		end

		zeta = scale * Ab{iblock}' * resid(:); % A' * W * (y - A*x)
		g = (rho * zeta + g) / (rho+1);

		% continuation
		rho = pi / (ii+1) * sqrt(1 - (pi / (2*ii+2))^2);
		ii = ii + 1;
	end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = single(x);
	end
	info(iter,:) = arg.userfun(x, arg.userarg{:});
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, varargin)
chat = evalin('caller', 'arg.chat');
if chat
%	x = evalin('caller', 'x');
	printm('minmax(x) = %g %g', min(x), max(x))
end
out = cpu('etoc');


% ir_pwls_os_lalm_test()
% test with small quadratic case 1/2 |y - A x|_2^2 + 1/2 |C x|_2^2
% see also ir_ct_fan_beam_sqs_vs_lalm.m
function ir_pwls_os_lalm_test
rng(0)
nd = 99; np = 38;
A = randn(nd,np);
A = abs(A);
R = Reg1(ones(np,1), 'type_penal', 'mat', 'type_diff', 'sparse', 'beta', 2^0);
C = R.C;
C = full(C);
if 0 % try to make it so unit step size is reasonable?
	D = abs(A)' * (abs(A) * ones(np,1)) + abs(C)' * (abs(C) * ones(np,1));
	A = A * diag(sqrt(1 ./ D));
	C = C * diag(sqrt(1 ./ D));
end
H = A' * A + C' * C;
pr cond(H)

xtrue = ones(np,1); xtrue(10:20) = 2;
yi = A * xtrue; yi = yi + 0.01 * max(yi(:)) * randn(nd,1);
xhat = H \ (A' * yi);
% nrms(xhat, xtrue)

x0 = zeros(np,1);
niter = 200;

printm 'Nesterov' % todo: replace with OGM
gradfun = @(x) A' * (A * x - yi) + C' * (C * x);
step = abs(A') * (abs(A) * ones(np,1)) + abs(C)' * (abs(C) * ones(np,1));
step = 1 / max(step); % todo: use D
xf = ir_fista(x0, 'gradfun', gradfun, 'step', step, 'niter', niter, ...
	'isave', 'all');

printm 'LALM'
xl = ir_pwls_os_lalm(x0, Gmatrix(A), yi, R, 'niter', niter, 'isave', 'all');

printm 'CG'
xc = qpwls_pcg1(x0, A, 1, yi, R.C, 'niter', niter, 'isave', 'all');

if im
	im plc 1 2
	im subplot 1
	plot([xtrue xhat xl(:,end) xc(:,end) xf(:,end)])
	im subplot 2
	xlabel 'iteration', ylabelf 'NRMSD $x^n$ vs $\hat{x}$'
	plot(	0:niter, nrms(xl, xhat), 'r-', ...
		0:niter, nrms(xf, xhat), 'g-', ...
		0:niter, nrms(xc, xhat), 'c-')
	legend('LALM', 'FGM', 'CG')
end
