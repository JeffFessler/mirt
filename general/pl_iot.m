 function [xs, info] = pl_iot(x, Ab, data, R, varargin)
%function [xs, info] = pl_iot(x, Ab, data, R, [options])
%|
%| Generic penalized-likelihood minimization,
%| for arbitrary negative log-likelihood with convex non-quadratic penalty,
%| via incremental optimization transfer using separable quadratic surrogates.
%|
%| cost(x) = sum_i h(data_i; [Ax]_i) + R(x),
%|
%| in
%|	x	[np 1]		initial estimate
%|	Ab	[nd np]		system matrix (see Gblock)
%|	data	{cell}		whatever data is needed for the likelihood
%|				data{1} must have size "[nb,na]", so for 3D
%|				one must use reshaper(yi, '2d').
%|	R			penalty object (see Robject.m)
%|
%| options
%|	dercurv	{function_handle} function returning derivatives and curvatures
%|				of negative log-likeihood via the call:
%|		[deriv curv] = dercurv(data, Ab{m}*x, curvtype, iblock, nblock)
%|				or choose 'trl' or 'wls' (default)
%|	niter 	#		# total iterations (default: 1+0, pure OS!)
%|	os	#		how many "warmup" OS iterations (default: 1)
%|	riter	#		# of penalty subiterations (default: 3).
%|	curvtype ''		curvature type:
%|					'pc' : precomputed (default)
%|					'oc' : optimal (ensures monotonic)
%|	pixmax	(value)		upper bound on pixel values (default: infinity)
%|	pixmin	(value)		lower bound on pixel values (default: 0)
%|	isave	[]		list of iterations to archive
%|				(default: [] 'last)
%|	gi	[nd]		precomputed Ab * 1 factors.
%|	chat	(value)		verbosity (default: 0)
%|	userfun			user defined function handle (see default below)
%|
%| out
%|	xs	[np nsave]	estimates each (saved) iteration
%|	info	[niter+1 ?]	userfun output.  default is cpu times
%|
%| Copyright 2005-3-8, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('pwls_example')
return
end
if nargin < 4, help(mfilename), error(mfilename), end

cpu etic

Ab = block_ob(Ab, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Ab, 'n');
starts = subset_start(nblock);

% defaults
arg.dercurv = 'wls';
arg.curvtype = 'pc';
arg.riter = 3;
arg.os = 1;
arg.niter = arg.os + 0; % default is pure OS, no IOT
arg.isave = [];
arg.pixmin = 0;
arg.pixmax = inf;
arg.gi = [];
arg.update_even_if_denom_0 = true;
arg.chat = false;
arg.userfun = @userfun_default;

% options
arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

if streq(arg.dercurv, 'wls')
	arg.dercurv = @wls_dercurv;
elseif streq(arg.dercurv, 'trl')
	arg.dercurv = @trl_dercurv;
end

if ~isa(arg.dercurv, 'function_handle')
	error 'dercurv not function handle?'
end

ticker(mfilename, 1, arg.niter)

[nb na] = size(data{1});

% g_i = sum_j g_ij.  caution: this requires real g_ij and g_ij >= 0
if isempty(arg.gi)
	arg.gi = sum(Ab');
	arg.gi = max(reale(arg.gi), 0); % trick: just in case...
	arg.gi = reshape(arg.gi, [nb na]);
end

% precompute likelihood-term denominator if needed
if streq(arg.curvtype, 'pc')
	[dummy curvi] = feval(arg.dercurv, data, 0, 'pc');
	ldenom = Ab' * col(arg.gi .* curvi); % one denominator shared by all subsets
else
	ldenoms = zeros(np, nblock);
end

np = length(x);
x = min(x, arg.pixmax);
x = max(x, arg.pixmin);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(arg.niter, ?); % trick: do not initialize since size may change

%
% precompute gradient-related state vectors, usually by OS-SPS warmup iterations
%
vvm = zeros(np,nblock);

if arg.os == 0
	for iblock=1:nblock
%		vvm(:,iblock) = Ab{ii}
		error 'not done'
	end
end


%
% SPS-OS iterations to "warm up" (often 1 or 2 iterations suffices)
%
for iter = 1:arg.os
	for iset = 1:nblock
		ticker([mfilename ': os'], [iter iset], [arg.os nblock])

		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Ab{iblock} * x;	% l=A*x "line integrals"
		li = reshape(li, nb, length(ia));
		[dhi curvi] = feval(arg.dercurv, data, li, ...
			arg.curvtype, iblock, nblock);

		if ~streq(arg.curvtype, 'pc') % on-the-fly curvatures
			ldenoms(:,iblock) = ...
				Ab{iblock}' * col(arg.gi(:,ia) .* curvi);
			ldenom = nblock * ldenoms(:,iblock);
		end

		lnum = ldenom .* x - nblock * (Ab{iblock}' * dhi(:));

		if iter == arg.os % save last gradient-related state vectors
			vvm(:,iblock) = lnum / nblock;
		end

		x = inner_update(x, lnum, ldenom, R, arg);
	end
	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
	info(1+iter,:) = feval(arg.userfun);
	if arg.chat, printf('Range %g %g', min(x), max(x)), end
end

%
% At this point we have initialized vv_m and ldenoms
% and initialized x based on the last subset.
% It is often logical to update x using all subsets now.
% This is "almost free" since no new likelihood gradients are used.
%
vv = sum(vvm,2);
if streq(arg.curvtype, 'pc')
	ldenoms = ldenom / nblock;
else
	ldenom = sum(ldenoms,2);
end
if 1
	x = inner_update(x, vv, ldenom, R, arg);
	if any(arg.isave == -1)
		xs(:, arg.isave == -1) = x;
	end
end


%
% IOT iterations
%
for iter = (1+arg.os):arg.niter
	for iset = 1:nblock
		ticker([mfilename ': iot'], [iter iset], [arg.niter nblock])

		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Ab{iblock} * x;	% l=A*x "line integrals"
		li = reshape(li, nb, length(ia));
		[dhi curvi] = feval(arg.dercurv, data, li, ...
			arg.curvtype, iblock, nblock);

		if ~streq(arg.curvtype, 'pc')
			ldenom = ldenom - ldenoms(:,iblock);
			ldenoms(:,iblock) = ...
				Ab{iblock}' * col(arg.gi(:,ia) .* curvi);
			ldenom = ldenom + ldenoms(:,iblock);
		end


		vv = vv - vvm(:,iblock);
		if streq(arg.curvtype, 'pc')
			vvm(:,iblock) = ldenoms .* x - (Ab{iblock}' * dhi(:));
		else
			vvm(:,iblock) = ldenoms(:,iblock) .* x - (Ab{iblock}' * dhi(:));
		end
		vv = vv + vvm(:,iblock);

		x = inner_update(x, vv, ldenom, R, arg);
	end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
	info(1+iter,:) = feval(arg.userfun);
	if arg.chat, printf('Range %g %g', min(x), max(x)), end
end


%
% inner_update()
% given SPS for likelihood, update x using regularization subiterations
%
function x = inner_update(x, lnum, ldenom, R, arg)

for ii=1:arg.riter
	if isempty(R)
		num = lnum;
		rdenom = 0;
	else
		rdenom = R.denom(R, x);
		num = lnum - R.cgrad(R, x) + rdenom .* x;
	end

	% update
	if arg.update_even_if_denom_0
		x = num ./ (ldenom + rdenom);
	else
		old = x(ldenom == 0);
		x = num ./ (ldenom + rdenom);
		x(ldenom == 0) = old;
	end
	x = max(x, arg.pixmin);	% lower bound
	x = min(x, arg.pixmax);	% upper bound
end


% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default
%x = evalin('caller', 'x');
out = [cpu('etoc')];
