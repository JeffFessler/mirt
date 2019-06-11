 function [xs, Ps] = qpwls_qn(x, G, W, yy, C, P, niter)
%function [xs, Ps] = qpwls_qn(x, G, W, yy, C, P, niter)
%
% quadratic penalized weighted least squares (QPWLS) via
% (preconditioned) quasi-Newton (QN) algorithm
% cost(x) = (y-Gx)'W(y-Gx) / 2  + x'C'Cx / 2
% in
%	x	[np,1]		initial estimate
%	G	[nn,np]		system matrix
%	W	[nn,nn]		data weighting matrix
%	yy	[nn,1]		noisy data
%	C	[nc,np]		penalty 'derivatives' (R = \Half C'*C)
%	P	[np,np]		preconditioner (matrix or object)
%	niter			# total iterations
% out
%	xs	[np,niter]	estimates each iteration
%	ws	[np,niter]	P rank-one update terms each iteration
%
% Copyright July 2000, Jeff Fessler, The University of Michigan

if nargin < 3 || nargin > 7, help(mfilename), error(mfilename), end
np = ncol(G);

if ~isvar('C')		| isempty(C),		C = 0;			end
if ~isvar('P')		| isempty(P),		P = speye(np);		end
if ~isvar('niter')	| isempty(niter),	niter = 2;		end
if ~isvar('x')		| isempty(x),		x = zeros(np,1);	end

yy = yy(:);
xs = zeros(np, niter);
xs(:,1) = x;

ws = zeros(np, niter);
bs = zeros(1, niter);

%
% initialize projections
%
Gx = G * x;
Cx = C * x;

grad_next = G' * (W * (yy - Gx)) - C' * Cx;

%
% iterate
%
for iter = 2:niter

	grad = grad_next;

	%
	% compute ddir = P * grad, which is -s
	%
	ddir = Pmul(grad, P, bs(1:(iter-2)), ws(:,1:(iter-2)));
	x = x + ddir;

	% check if descent direction
	if ddir' * grad < 0
		warning('wrong direction')
%		keyboard
	end

	%
	% next (negative) gradient
	%
	Gx = G * x;
	Cx = C * x;
	grad_next = G' * (W * (yy - Gx)) - C' * Cx;

	q = grad_next - grad;
	w = ddir + Pmul(q, P, bs(1:(iter-2)), ws(:,1:(iter-2)));

	tmp = w' * q;
	if tmp == 0, error zero, end
	bs(iter-1) = 1 / tmp;
	ws(:,iter-1) = w;

	%
	% step size in search direction
	%
	if 0
		Gdir = G * ddir;
		Cdir = C * ddir;

		step = (ddir' * grad) / (Gdir'*(W*Gdir) + Cdir'*Cdir);
		steps(iter,1) = step;
		if step < 0
			warning('downhill?')
			keyboard
		end

		%
		% update
		%
		x	= x + step * ddir;
		Gx	= Gx  + step * Gdir;
		Cx	= Cx  + step * Cdir;
	end

	xs(:,iter) = x;
end


%
% multiply v = P_n * u where P_n = P_0 + \sum_i b_i w_i w_i'
%
function v = Pmul(u, P, bs, ws)
v = P * u;
bs = bs .* (u' * ws);
v = v - ws * bs';
