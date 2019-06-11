 function xs = eml_osem_mismatch ...
	(x, Gf, Gback, yi, ci, ri, niter, pixmax, precon, relax0)
% E-ML-OSEM algorithm for image reconstruction from Poisson emission data
% (ordered subsets expectation maximization)
% model: Y_i ~ Poisson(c_i [G x]_i + r_i)
% in
%	x	[np,1]		initial estimate
%	Gf	[nd,np]		Gblock forward projector
%	Gback	[nd,np]		Gblock back projector (i.e., mismatched)
%	yi,ci,ri [nb,na]	see em_fbp.m (for model too)
%	niter			# iterations
%	pixmax			upper constraint for pixel values
%	precon [np,{1|nblock}]	preconditioners
%	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
% out
%	xs [np,niter]	updated image vectors each iteration
%
% Copyright 2002-3-14, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

Gf = block_ob(Gf, 'ensure'); % make it a block object (if not already)
nblock = Gf.nblock;
starts = subset_start(nblock);

if ~isa(Gback, 'Gblock') || nblock ~= Gback.nblock
	error 'Gback must be Gblock with same number of blocks'
end

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('niter')	| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	| isempty(pixmax),	pixmax = inf;	end
if ~isvar('chat')	| isempty(chat),	chat = false; end

if ~isvar('relax0')	| isempty(relax0),	relax0 = 1;	end
if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

eml_check(yi, ci, ri, 'os', nblock);

[nb, na] = size(yi);

flag_classic = false;
if ~isvar('precon') || isempty(precon)
	warning 'Using fast preconditioner rather than classic OSEM'
	precon = nblock ./ (Gf' * ci(:));
elseif streq(precon, 'classic')
	flag_classic = true;
else	% user-specified precon
	if ncol(precon) ~= nblock && ncol(precon) ~= 1
		error 'precon columns must be 1 or nblock'
	end
end


%
% loop over iterations
%
xs = zeros(numel(x), niter);
x = max(x,0);
x = min(x,pixmax);
xs(:,1) = x(:);

for iter = 2:niter
	if chat, printf('E-ML-OS-EM iteration %d', iter-1), end

	relax = relax0 / (1 + relax_rate * (iter-2));

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		% predicted measurements
		yp = reshape(Gf{iblock} * x, nb, length(ia));
		yp = ci(:,ia) .* yp + ri(:, ia);
		yp(yp == 0) = inf;	% avoids /0 error

		if flag_classic
			Asum = Gback{iblock}' * col(ci(:,ia));
			pre = ones(size(Asum));
			pre(Asum > 0) = 1 ./ Asum(Asum > 0);
		else
			pre = precon(:,min(iset,ncol(precon)));	% 1 or iset
		end

		dhi = ci(:,ia) .* (yi(:,ia) ./ yp - 1);
		grad = Gback{iblock}' * dhi(:);

		x = x + relax * (x .* pre) .* grad;
		x = max(x,0);
		x = min(x,pixmax);
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,iter) = x;
end
