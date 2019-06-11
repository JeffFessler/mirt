 function xs = eml_inc_em(x, Gb, yi, ci, ri, varargin)
%|function xs = eml_inc_em(x, Gb, yi, ci, ri, varargin)
%| E-ML-INC-EM algorithm for image reconstruction from Poisson emission data
%| (incremental EM algorithm)
%| model: Y_i ~ Poisson(c_i [G x]_i + r_i)
%| in
%|	x	[np,1]		initial estimate
%|	Gb	[nd,np]		Gblock object (see eml_osem_test.m)
%|	yi,ci,ri [nb,na]	see em_fbp.m (for model too)
%|
%| options
%|	'niter'	(#)		# of iterations
%|	'isave'	[]		list of iterations to archive
%|				(default: [] 'last)
%|	chat	(0 | 1)		verbosity
%|	pixmax	(value)		upper constraint for pixel values
%|	hds	(value)		'1' or '3' hidden data space
%|	os	(#)		how many "warmup" osem iterations
%|
%| out
%|	xs [np nsave]		estimates each (saved) iteration
%|
%| Copyright 2004-3-20, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

% options
arg.chat = 0;
arg.pixmax = inf;
arg.hds = 1;
arg.niter = 1;
arg.isave = 'last';
arg.os = 0;
arg.userfun = [];

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

Gb = block_ob(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

ticker(mfilename, 1, arg.niter)

% precompute hidden data space factor
gam = eml_hds(Gb, ci, ri, arg.hds);
if arg.chat, printf('hds = %g', gam), end

% system (block) sensitivities
np = length(x);
asum = zeros(np, 1);
precon = zeros(np, nblock); % classic OSEM preconditioner
for iset=1:nblock
	ticker([mfilename ' : precon'], iset, nblock)
	istart = starts(iset);
	ia = istart:nblock:na;
	asum_m = Gb{istart}' * col(ci(:,ia));
	asum = asum + asum_m;
	asum_m(asum_m == 0) = Inf; % avoid divide by 0
	precon(:, iset) = 1 ./ asum_m;
end, clear asum_m
% asum = Gb' * ci(:);

x = max(x,0);
x = min(x,arg.pixmax);
xs = zeros(numel(x), length(arg.isave));
if any(arg.isave == 0)
	xs(:,find(arg.isave == 0)) = x;
end


%
% precompute sufficient statistics,
% possibly using osem warmup iterations
%
ff = zeros(np,nblock);
if arg.os == 0
	for iblock=1:nblock
		ff(:,iblock) = compute_fm(x, Gb, yi, ci, ri, gam, iblock);
	end

else
	for iter = 1:arg.os
		for iset = 1:nblock
			ticker([mfilename ': osem'], [iter iset], [arg.os nblock])
			iblock = starts(iset);
			ia = iblock:nblock:na;

			tmp = compute_fm(x, Gb, yi, ci, ri, gam, iblock);

			pre = precon(:,min(iset,ncol(precon)));	% 1 or iset
			x = max(pre .* tmp - gam, 0);

			if iter == arg.os
				ff(:,iblock) = tmp;
			end
		end

		if any(arg.isave == iter)
			xs(:, arg.isave == iter) = x;
		end
	end
end
fsum = sum(ff, 2);


%
% loop over incremental EM iterations
%

for iter = (2+arg.os):arg.niter

	%
	% loop over subsets
	%
	for iset=1:nblock
		ticker(mfilename, [iter iset], [arg.niter nblock])
		iblock = starts(iset);

		x = fsum ./ asum - gam;
		x = max(x,0);
		x = min(x,arg.pixmax);

		fsum = fsum - ff(:,iblock);
		ff(:,iblock) = compute_fm(x, Gb, yi, ci, ri, gam, iblock);
		fsum = fsum + ff(:,iblock);
	end

	if arg.chat, printf('Range %g %g', min(x), max(x)), end
	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
end


%
% compute_fm()
% compute sufficient statistics
%
function fm = compute_fm(x, Gb, yi, ci, ri, gam, iblock)
eterm = eml_eterm(x, Gb, yi, ci, ri, iblock);
fm = (x + gam) .* eterm;
