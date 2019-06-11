 function [xs, info] = ir_diag_majorize_admm1(H, varargin)
%function [xs, info] = ir_diag_majorize_admm1(H, varargin)
%|
%| Use ADMM to find a diagonal matrix D = diag(x) of minimum Frobenious norm
%| that majorizes a given square matrix H, i.e., min_D || D || s.t. H <= D.
%| Assumes that H is positive semidefinite (so thus D is as well).
%|
%| Solved by min_{D, S >=0} || D ||_F s.t. S = D - H
%|
%| in
%|	H	[n n]		matrix to be majorized.
%|
%| options
%|	init	[n]		initial guess of diagonal (default: H * 1)
%|	mu	[1]		AL penalty parameter (default: 1)
%|
%|	niter			# total iterations (default: 1)
%|					(max # if tol used)
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun	@		user defined function handle (see default below)
%|					taking arguments (x, iter, userarg{:})
%|	userarg {}		user arguments to userfun (default {})
%|	stop_diff_tol		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 0
%|	stop_diff_norm		use norm(., type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	d	[n]		diagonal of majorizing diagonal matrix
%|
%| Copyright 2013-11-02, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage(), end
if nargin == 1 && streq(H, 'test'), ir_diag_majorize_admm1_test, return, end

% defaults
arg.init = '';
arg.mu = 1;

arg.niter = 1;
arg.isave = [];
arg.vnorm = 2; % default is frobeenius norm squared, i.e., |vec(D)|_2
arg.userfun = @userfun_default;
arg.userarg = {};
arg.stop_diff_tol = 0;
arg.stop_diff_norm = 2;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.init)
	arg.init = sum(H, 2); % H * 1
end

mu = arg.mu;

x = arg.init;
np = numel(x);

if ~isequal(size(H), [np np]) || ~isequal(size(x), [np 1])
	fail 'size of H and initial x must match'
end

arg.isave = iter_saver(arg.isave, arg.niter);
if arg.stop_diff_tol
	norm_diff = @(x) norm(x(:), arg.stop_diff_norm);
end

cpu etic

xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end


E = 0;
S = diag(x) - H;

softshrink = @(b, reg) sign(b) .* max(abs(b) - reg, 0);

% iterate
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	xp = x; % previous

	% update primal
	switch arg.vnorm
	case 2 
		x = mu / (1+mu) * diag(S + H - E);
	case 1
		x = diag(S + H - E);
		x = softshrink(x, 1/mu); % 1/2 | S - (D-H) - E |^2 + 1/mu |D|_1
	otherwise
		fail('bad arg.vnorm %g', arg.vnorm)
	end

	if any(x < 0)
		minmax(x)
		x = max(x, 0); % keep D positive semidefinite, i.e., D >= 0
	end
	D = diag(x);

	% update auxiliary
	tmp = D - H + E;
	[V, lam] = eig(tmp);
	lam = diag(lam);
	lam = reale(lam);
	lam = max(lam, 0); % hard thresholding of eigenvalues
	S = V * diag(lam) * V';

	% update multiplier
	E = E - (S + H - D);

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x(:);
	end
	info(iter,:) = arg.userfun(x, iter, arg.userarg{:});
%	pr 'info(iter,1)'

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if iter > 1 && arg.stop_diff_tol
		ratio = norm_diff(x - xp) / norm_diff(x);
		if ratio < arg.stop_diff_tol
			if 1 || arg.chat
				printm('stop at iteration %d, diff %g < %g', ...
					iter, ratio, arg.stop_diff_tol)
			end
			if isequal(arg.isave, arg.niter) % saving only last?
				xs = x(:); % save the 'final' iterate
			else % saving many iterates?
				xs(:, arg.isave > iter) = []; % clear out unused
			end
		break
		end
	end
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default(x, iter, varargin)
%pr minmax(x)
out = [cpu('etoc')];


% ir_diag_majorize_enforce()
% find (smallest) alpha such that H <= diag(x) + alpha I
% and return x + alpha, to enforce majorization
function x = ir_diag_majorize_enforce(x, H)
lam = eig(diag(x) - H);
x = x - min(min(lam), 0);


function ir_diag_majorize_admm1_test
n = 2^3;
A = Gblur(true(n,1), 'psf', [0 1 -1]'); % minor reduction
%A = randn(n+3,n); % this case gives dramatic reduction!
H = full(A' * A);
x0 = abs(A)' * (abs(A) * ones(n,1)); % per SQS theory

% vnorm = 1; name = '1-'; % no benefit in this example
vnorm = 2; name = 'Frobenius';
niter = 100;
xs = ir_diag_majorize_admm1(H, 'init', x0+0.0, 'mu', 2^8, ...
	'vnorm', vnorm, ...
	'niter', niter, 'isave', 'all');

frob = sum(abs(xs).^vnorm, 1);
frob = frob / norm(x0,vnorm)^vnorm; % normalize to [0,1]

lam = zeros(n, niter+1);
for ii=1:niter+1
	x = xs(:,ii);
	lam(:,ii) = eig(diag(x) - H);
	xe(:,ii) = ir_diag_majorize_enforce(xs(:,ii), H); % enforced
end
equivs(xe(:,1), xs(:,1))

frobe = sum(abs(xe).^vnorm, 1);
frobe = frobe / norm(x0,vnorm)^vnorm;

im plc 2 2

im subplot 1
plot(0:niter, frob, '.-')
xlabel 'iteration'
titlef('Relative %s norm', name)
axis([0 niter 0.98 1.01])
ytick([0.98 1 1.01])

im subplot 2
plot(0:niter, min(lam), '.-')
titlef('Minimum eigenvalue of $D-H$')
axis([0 niter -0.04 0.03])
ytick([-0.04 0 0.03])

im subplot 3
%plot(0:niter, max(lam), '-o')
%plot(0:niter, sum(xs), '-o')
plot(0:niter, frobe, '.-')
titlef('Relative %s norm of $D+\alpha I$', name)
ir_text(.6*niter, 1.002, '$\alpha = -\min(\lambda_{min}(D-H),0)$', ...
	'horizontalalignment', 'center')
xlabelf('iteration')
ymin = min(0.994, min(frobe));
ymax = max(1.008, max(frobe));
axis([0 niter ymin ymax])
ytick([ymin 1 ymax])

im subplot 4
plot(1:n, x0, '-o', 1:n, xe(:,end), '-x')
ir_legend({'initial', 'final'}, 'location', 'south')
titlef('Diagonal values')
xlabelf '$j$'
ylabelf '$D_{jj}$'
axis([1 n 1.9 4.1])
ytick([2 3 4])

if any(eig(diag(xe(:,end)) - H) < 0), fail 'bug', end

% ir_savefig cw ir_diag_majorize_admm1
