  function xs = eql_os_emdp(x, Gb, yi, ci, ri, R, varargin)
%|function xs = eql_os_emdp(x, Gb, yi, ci, ri, R, [options])
%|
%| Ordered-subsets version of Alvaro De Pierro's modified EM algorithm
%| for quadratically penalized log-likelihood emission image reconstruction
%| from Poisson data.
%| Based on Mar. 1995 IEEE T-MI paper: A modified expectation maximization...
%|
%| see eql_sps_os.m for arguments
%|
%| options
%|	niter	(#)		# of iterations (default: 1)
%|	isave	[array]		which iterations to archive (default: 'all')
%|	relax0			initial relaxation (default: 1)
%|				and (optionally) decay rate (default: 0)
%|	pixmax	(value)		upper bound on pixel values (default: infinity)
%|	hds	(1 or 3)	which hidden data space (default: 1)
%|	chat	(value)		verbosity (default: 0)
%
%| Copyright 2002-3-14, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.niter = 1;
arg.pixmax = inf;
arg.chat = false;
arg.hds = 1;
arg.relax0 = 1;
arg.isave = 'all';
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

if length(arg.relax0) == 1
	arg.relax_rate = 0;
elseif length(arg.relax0) == 2
	arg.relax_rate = arg.relax0(2);
	arg.relax0 = arg.relax0(1);
else
	error relax
end

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

ticker(mfilename, 0, arg.niter)

% precompute hidden data space factor
gam = eml_hds(Gb, ci, ri, arg.hds);
if arg.chat, printf('hds = %g', gam), end

%
% precompute the partial matrix sums (ala classic OSEM)
%
for iset=1:nblock
	ticker([mfilename ' : precon'], iset, nblock)
	istart = starts(iset);
	ia = istart:nblock:na;
	aj_sets(:,iset) = Gb{istart}' * col(ci(:,ia));
end
if arg.chat
printf('%d of %d zeros in aj_sets', sum(aj_sets(:) == 0), numel(aj_sets))
end

% for unpenalized case
if isempty(R)
	pgrad = 0;
	Rdenom = 0;
	aj_sets(aj_sets == 0) = Inf; % avoid divide by 0 later
end

%
% loop over iterations
%
xs = zeros(numel(x), length(arg.isave));
if any(x <= 0), error 'need x > 0', end
x = min(x, arg.pixmax);
if any(arg.isave == 0)
	xs(:,find(arg.isave == 0)) = x;
end

for iter = 1:arg.niter
	relax = arg.relax0 / (1 + arg.relax_rate * (iter-1));

	%
	% loop over subsets
	%
	for iset=1:nblock
		ticker(mfilename, [iter iset], [arg.niter nblock])

		iblock = starts(iset);
		eterm = eml_eterm(x, Gb, yi, ci, ri, iblock);

		if ~isempty(R)
			pgrad = R.cgrad(R, x);
			Rdenom = R.denom(R, x);
		end

		% trick: scale by nblock, to maintain resolution, per Aspire
		aj = nblock * aj_sets(:,iset);
		ej = nblock * eterm;
		gam_os = gam / nblock;	% trick: subset gamma factor

		xold = x;
		x = eql_root(Rdenom, ...
			(aj + pgrad - Rdenom .* (x+gam_os)) / 2, ...
			(x+gam_os) .* ej) - gam_os;
		x = (1-relax) * xold + relax * x;
		x = max(x, 0);
		x = min(x, arg.pixmax);
	end

	if arg.chat, printf('Range %g %g', min(x), max(x)), end

	if any(arg.isave == iter)
		xs(:,find(arg.isave == iter)) = x;
	end
end
