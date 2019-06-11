 function xs = epl_inc(x, Gb, yi, ci, ri, R, varargin)
%function xs = epl_inc(x, Gb, yi, ci, ri, R, [options])
%| Incremental optimization transfer version of EM algorithm for
%| penalized log-likelihood emission image reconstruction from Poisson data.
%| see eql_sps_os.m for arguments
%| in
%|	Gb			cell array or block object; see block_op()
%|	yi,ci,ri		cell arrays of vectors; see block_op()
%|
%| options
%|	niter	(#)		# of iterations (default: 1+1)
%|	isave	[]		which iterations to archive (default: 'all')
%|	relax0			initial relaxation (default: 1)
%|				and (optionally) decay rate (default: 0)
%|	starts	[nblock]	order of view starts 
%|	pixmax	(value)		upper bound on pixel values (default: infinity)
%|	hds	(1 or 3)	which hidden data space (default: 1)
%|	os	(#)		how many "warmup" OSDP iterations (default: 1)
%|	chat	(value)		verbosity (default: 0)
%|
%| Copyright 2005-2-17, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

Gb = block_op(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_op(Gb, 'n');

% defaults
arg.niter = 2;
arg.isave = 'all';
arg.pixmax = inf;
arg.chat = false;
arg.hds = 1;
arg.os = 1;
arg.starts = subset_start(nblock); % fix: this is a kludge!
arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);

yi = block_op(Gb, 'ensure_block_data', yi);

if ~isvar('ci') || isempty(ci)
	ci = num2cell(ones(1,nblock));
else
	ci = block_op(Gb, 'ensure_block_data', ci);
end

if ~isvar('ri') || isempty(ri)
	ri = num2cell(zeros(1,nblock));
else
	ri = block_op(Gb, 'ensure_block_data', ri);
end

eml_check(yi, ci, ri, 'os', nblock);

if ~isvar('R'), R = []; end

ticker(mfilename, 1, arg.niter)

% precompute hidden data space factor
gam = eml_hds(Gb, ci, ri, arg.hds);
if arg.chat, printf('hds = %g', gam), end


% system (block) sensitivities
x = x(:);
np = length(x);
asum = zeros(np, 1);
aj_sets = cell(1,nblock);
for iblock=1:nblock
	ticker([mfilename ' : aj'], iblock, nblock)
	istart = arg.starts(iblock);
	aj_sets{istart} = Gb{istart}' * ci{istart};
	asum = asum + aj_sets{istart};
end

if arg.chat
	tmp = cat(1, aj_sets{:});
	printf('%d of %d zeros in aj_sets', sum(tmp == 0), numel(tmp))
	clear tmp
end

% for unpenalized case
if isempty(R)
	pgrad = 0;
	Rdenom = zeros(size(x));
	for iblock=1:nblock
		aj_sets{iblock}(aj_sets{iblock} == 0) = Inf; % avoid / 0 later
	end
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
		ff(:,iblock) = (x+gam) .* eml_eterm(x, Gb{iblock}, ...
				yi{iblock}, ci{iblock}, ri{iblock});
	end

else
	for iter = 1:arg.os
		for iblock = 1:nblock
			ticker([mfilename ': os'], [iter iblock], [arg.os nblock])

			istart = arg.starts(iblock);
			eterm = eml_eterm(x, Gb{istart}, ...
					yi{istart}, ci{istart}, ri{istart});

			if iter == arg.os
				ff(:,iblock) = (x + gam) .* eterm;
			end

			if ~isempty(R)
				pgrad = R.cgrad(R, x);
				Rdenom = R.denom(R, x);
			end

			% trick: scale by nblock, for resolution, per Aspire
			aj = nblock * aj_sets{istart};
			ej = nblock * eterm;
			gam_os = gam / nblock; % trick: subset gamma factor

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
	for iblock=1:nblock
		ticker(mfilename, [iter iblock], [arg.niter nblock])

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

		istart = arg.starts(iblock);
		eterm = eml_eterm(x, Gb{istart}, ...
				yi{istart}, ci{istart}, ri{istart});

		fsum = fsum - ff(:,iblock);
		ff(:,iblock) = (x + gam) .* eterm;
		fsum = fsum + ff(:,iblock);
	end

	if arg.chat, printf('Range %g %g', min(x), max(x)), end

	if any(arg.isave == iter)
		xs(:,find(arg.isave == iter)) = x;
	end
end
