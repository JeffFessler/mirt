 function [xs, info] = qpwls_pcg(x, G, W, yi, nder1, C, M, niter, mask, chat)
%function [xs, info] = qpwls_pcg(x, G, W, yi, nder1, C, M, niter, mask, chat)
%
% quadratic penalized weighted least squares via
% preconditioned conjugate gradients algorithm
% cost:	Psi(x) = (y-Gx)'W(y-Gx)/2 - n'(y-Gx) + x'C'Cx/2
% in
%	x	[np,1] or [nx,ny]	initial estimate
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix
%	yi	[nd,1]		noisy data
%	nder1	[nd,1]		linear term (obsolete: just use zero!)
%	C	[nc,np]		penalty 'derivatives' (R = \Half C'*C)
%	M	[np,np]		preconditioner (or object) (or 1)
%				use 'circ0' for center-based circulant precon!
%	niter			# total iterations
%	mask	[nx,ny]		caution: used only to initialize and display!
% out
%	xs	[np,niter]	estimates each iteration
%	info	[niter, 4]	gamma, step, time, flops
%
% Copyright Jan 1998, Jeff Fessler, The University of Michigan

if nargin < 8 || nargin > 10, help(mfilename), error(mfilename), end

if ~isvar('chat') || isempty(chat), chat = 0; end

if ~isreal(yi) && ~isempty(M)
	persistent warned
	if isempty(warned), warned = 0; end
	if ~warned
		warning 'not 100% sure about the complex preconditioned case'
		warned = 1;
	end
end

%
% if requested, build circulant preconditioner based on center pixel
%
if ischar(M) && streq(M, 'circ0')
	if ~isvar('mask') || isempty('mask')
		error 'mask required for circulant precon'
	end
	M = qpwls_precon(M, {G, W}, C, mask);
end

np = ncol(G);
if np ~= 1
	if numel(x) ~= np && isvar('mask')
		x = x(mask(:));
	end
else
	x = x(:);
	np = length(x);
end

yi = yi(:);
xs = zeros(np, niter);
xs(:,1) = x;

if chat > 1
	im clf, im(121, embed(x,mask), 'xinit')
end

info = zeros(niter,4);

%
% initialize projections
%
Gx = G * x;
Cx = C * x;

tic
%flops0 = flops;

%
% iterate
%
for iter = 2:niter
	ticker(mfilename, iter, niter)

	%
	% (negative) gradient
	%
	ngrad = G' * (W * (yi-Gx) - nder1) - C' * Cx;

	%
	% preconditioned gradient
	%
	pregrad = M * ngrad;

	% todo: constrain pregrad to be real, if user wants x to be real

	%
	% direction
	%
	newinprod = ngrad' * pregrad;
	% fix: should i take the real part?
	newinprod = reale(newinprod, 'warn', 'inprod');
	if iter == 2
		ddir = pregrad;
		gamma = 0;
	else
	%	gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		if (oldinprod == 0)
			warning 'inprod=0.  going nowhere!'
			gamma = 0;
		else
			gamma = newinprod / oldinprod;	% Fletcher-Reeves
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if real(ddir' * ngrad) < 0
		warning 'wrong direction'
		keyboard
	end

	%
	% step size in search direction
	%
	Gdir = G * ddir;
	Cdir = C * ddir;

	denom = Gdir'*(W*Gdir) + Cdir'*Cdir;
	if denom == 0
		warning 'found exact solution???  step=0 now!?'
		step = 0;
	else
		denom = reale(denom, 'error', 'denom');
		step = (ddir' * ngrad) / denom;
%		step = reale(step, 'warn', 'step');
		step = real(step); % real step sizes seems only logical
	end
	if step < 0
		warning 'downhill?'
		keyboard
	end

	%
	% update
	%
	Gx	= Gx  + step * Gdir;
	Cx	= Cx  + step * Cdir;
	x	= x + step * ddir;
	xs(:,iter) = x;

	if chat > 1 && (iter < 10 || rem(iter,10)==0)
		im(121, embed(x,mask), 'x')
		title(sprintf('%d step=%g gamma=%g', iter, step, gamma))
		if isobject(M)
			ptype = M.type;
		else
			ptype = 'unknown';
		end
		im(122, embed(pregrad,mask), ['PreGrad ' ptype])
		drawnow
	end

	info(iter,1) = gamma;
	info(iter,2) = step;
	info(iter,3) = toc;		% accum. time
%	info(iter,4) = flops;		% accum. flops
end

%	info(:,4) = info(:,4) - flops0;
