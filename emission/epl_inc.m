 function xs = epl_inc(x, Gb, yi, ci, ri, R, varargin)
%function xs = epl_inc(x, Gb, yi, ci, ri, R, [options])
% Incremental optimization transfer version of EM algorithm for
% penalized log-likelihood emission image reconstruction from Poisson data.
% see eql_sps_os.m for arguments
%
% options
%	niter	(#)		# of iterations (default: 1+1)
%	isave	[]		which iterations to archive (default: 'all')
%	relax0			initial relaxation (default: 1)
%				and (optionally) decay rate (default: 0)
%	pixmax	(value)		upper bound on pixel values (default: infinity)
%	hds	(1 or 3)	which hidden data space (default: 1)
%	os	(#)		how many "warmup" OSDP iterations (default: 1)
%	chat	(value)		verbosity (default: 0)
%
% Copyright 2005-2-17, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.niter = 2;
arg.isave = 'all';
arg.pixmax = inf;
arg.chat = false;
arg.hds = 1;
arg.os = 1;
arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

Gb = block_ob(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi(:)));
end

if ~isvar('R'), R = []; end

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

ticker(mfilename, 1, arg.niter)

% precompute hidden data space factor
gam = eml_hds(Gb, ci, ri, arg.hds);
if arg.chat, printf('hds = %g', gam), end


% system (block) sensitivities
np = length(x);
asum = zeros(np, 1);
%precon = zeros(np, nblock); % classic OSEM preconditioner
for iset=1:nblock
	ticker([mfilename ' : aj'], iset, nblock)
	istart = starts(iset);
	ia = istart:nblock:na;
	aj_sets(:,iset) = Gb{istart}' * col(ci(:,ia));
end
asum = sum(aj_sets,2);
if arg.chat
printf('%d of %d zeros in aj_sets', sum(aj_sets(:) == 0), numel(aj_sets))
end

% for unpenalized case
if isempty(R)
	pgrad = 0;
	Rdenom = zeros(size(x));
	aj_sets(aj_sets == 0) = Inf; % avoid divide by 0 later
end

if any(x < 0), error 'need x >= 0', end
x = min(x,arg.pixmax);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:,find(arg.isave == 0)) = x;
end

%
% precompute sufficient statistics,
% probably using OSDP warmup iterations
%
ff = zeros(np,nblock);
if arg.os == 0
	'warning: not tested'
	for iblock=1:nblock
		ff(:,iblock) = (x+gam) .* eml_eterm(x, Gb, yi, ci, ri, iblock);
	end

else
	for iter = 1:arg.os
		for iset = 1:nblock
			ticker([mfilename ': os'], [iter iset], [arg.os nblock])

			iblock = starts(iset);
			eterm = eml_eterm(x, Gb, yi, ci, ri, iblock);

			if iter == arg.os
				ff(:,iblock) = (x + gam) .* eterm;
			end

			if ~isempty(R)
				pgrad = R.cgrad(R, x);
				Rdenom = R.denom(R, x);
			end

			% trick: scale by nblock, for resolution, per Aspire
			aj = nblock * aj_sets(:,iset);
			ej = nblock * eterm;
			gam_os = gam / nblock;  % trick: subset gamma factor

			x = eql_root(Rdenom, ...
				(aj + pgrad - Rdenom .* (x + gam_os)) / 2, ...
				(x + gam_os) .* ej) - gam_os;
			x = max(x, 0);
			x = min(x, arg.pixmax);
		end
		if any(arg.isave == iter)
			xs(:,find(arg.isave == iter)) = x;
		end
	end, clear aj ej gam_os aj_sets
end
fsum = sum(ff, 2);

%
% loop over iterations
%
for iter = (1+arg.os):arg.niter

	%
	% loop over subsets
	%
	for iset=1:nblock
		ticker(mfilename, [iter iset], [arg.niter nblock])

		if ~isempty(R)
			pgrad = R.cgrad(R, x);
			Rdenom = R.denom(R, x);
		end

		% fix: inner subiterations for non-quadratic penalties?
		x = eql_root(Rdenom, ...
			(asum + pgrad - Rdenom .* (x + gam)) / 2, ...
			fsum) - gam;
		x = max(x, 0);
		x = min(x, arg.pixmax);

		iblock = starts(iset);
		eterm = eml_eterm(x, Gb, yi, ci, ri, iblock);

		fsum = fsum - ff(:,iblock);
		ff(:,iblock) = (x + gam) .* eterm;
		fsum = fsum + ff(:,iblock);
	end

	if arg.chat, printf('Range %g %g', min(x), max(x)), end

	if any(arg.isave == iter)
		xs(:,find(arg.isave == iter)) = x;
	end
end
