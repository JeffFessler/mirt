 function [ws, steps] = samsonov(wo, niter, J, W, M)
%function [ws, steps] = samsonov(wo, niter, J, W, M)
%
% Samsonov's DCF method:
% quadratic (weighted) least squares (QWLS) via
% preconditioned steepest descent (PSD) algorithm,
% over non-negative wi's only 
% cost(w) = (u-Jw)'W(u-Jw) / 2 
% in
%	wo	[nd,1]		initial estimate
%	niter			# total iterations
%   J	[nd,nd]		B*WB, Fatrix approx of Fourier crosstalk matrix, K
%
%   optional inputs:
%	W	[nn,nn]		data weighting matrix, default is I
%	M	[nd,nd]		preconditioning (diagonal) matrix  
%                   M = diag(wJack) adds stability and is the default,
%                   where wJack are Jackson's density compensation weights
% out
%	ws	[nd,niter]	density compensation weights computed at each iteration
%	steps	[niter,1]	step size each iteration
%
% NOT YET TESTED for W ~= I, M ~= diag(wJack)
%
% Kim Khalsa, Nov 2005, modified from:
% qpwls_psd.m, Copyright Jun 2000, Jeff Fessler, The University of Michigan

   
if nargin < 3 | nargin > 6, help(mfilename), error(mfilename), end
nd = length(wo);

if ~isvar('niter') | isempty(niter), niter = 2; end
if ~isvar('W'), W = eye(nd); end
if ~isvar('u'), u = ones(size(wo)); end
if ~isvar('M'), M = diag(wo); end

ws = zeros(nd, niter);
ws(:,1) = wo;
w = wo;
z = zeros(size(w));

steps = zeros(niter,2);

%
% iterate
%
for ii=2:niter
	%
	% (negative) gradient
	%

    grad = J' * (W * (u - (J * w)));

	%
	% preconditioned gradient: search direction
	%
	if ~isempty(M)
		ddir = M * grad;
	else
		ddir = grad;
    end

 	Jdir = J * ddir;

	% check if descent direction
	if ddir' * grad < 0
		warning('wrong direction')
		keyboard
	end

	%
	% step size in search direction
	%
	step = (ddir' * grad) / (Jdir'*(W*Jdir));
	steps(ii,1) = step;
	if step < 0
		warning('downhill?')
		keyboard
	end

	%
	% update
	%
	w	= w + step * ddir;
    w_tmp = w;
    w   = max([real(w) z]')';    % wi's must be non-negative
	ws(:,ii) = w;
end