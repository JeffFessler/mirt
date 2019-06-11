  function xs = tml_os_mltr(x, Ab, yi, bi, ri, varargin)
%|function xs = tml_os_mltr(x, Ab, yi, bi, ri, varargin)
%|
%| T-ML-TR algorithm of Nuyts et al for transmission Poisson problem
%| (PMB Apr. 1998 doi 10.1088/0031-9155/43/4/003)
%|
%| in
%|	x		[np 1]		initial guess
%|	Ab		[nd np]		system matrix (or Gblock object)
%|	yi,bi,ri	[nb na]		see tr_fbp.m (for model too)
%|
%| option 
%|	'niter'		# of iterations (default: 1)
%|	'isave'		which iterations to save (default: 'last')
%|	'pixmin'	lower bound for pixel values (default 0)
%|	'pixmax'	upper bound for pixel values (default inf)
%|				(these can be scalar or arrays the size of x)
%|	'precon'	'classic' (default) or 'fast' or 'nuyts'
%|	'ai'		[nb na]	precomputed 1' * Ab
%|	'userfun'	function to be called after every update if provided
%|				can use: x = evalin('caller', 'x');
%|
%| out
%|	xs [np niter]	updated image vectors each (saved) iteration
%|
%| Copyright 2010-07-11, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi), class(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi), class(yi));
end

% options
arg.userfun = [];
arg.pixmin = 0;
arg.pixmax = inf;
arg.isave = 'last';
arg.niter = 1;
arg.step = 1;
arg.precon = 'classic'; % or 'fast'
arg.ai = []; % sum of each row

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

ticker(mfilename, 1, arg.niter)

trl_check(yi, bi, ri);

% sum of each row
if isempty(arg.ai)
	arg.ai = reshape(sum(Ab'), size(yi)); % [nd 1] a_i = sum_j a_ij
end

Ab = block_ob(Ab, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Ab, 'n');
starts = subset_start(nblock);

switch arg.precon % denominator
case 'nuyts'
	aimax = max(arg.ai(:));
case 'classic'
	if any(ri(:) ~= 0)
		warn 'untested for ri=0'
	end
case 'fast'
	ldenom = Ab' * col((1 - div0(ri, yi)) .* (yi - ri) .* arg.ai);
otherwise
	fail('bad precon: %s', arg.precon)
end

% initialize
np = length(x);
x = min(x, arg.pixmax);
x = max(x, arg.pixmin);
xs = zeros(np, length(arg.isave), class(yi));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

[nb na] = size(yi);

for iter=1:arg.niter % iterations
	for iset = 1:nblock % subsets
		ticker([mfilename ': os'], [iter iset], [arg.niter nblock])

		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Ab{iblock} * x; % l=A*x "line integrals"
		li = reshape(li, nb, length(ia));

		bel = bi(:,ia) .* exp(-li);
		yb = bel + ri(:,ia); % predicted measurement means 

		dothi = (1 - yi(:,ia) ./ yb) .* (-bel);

		% negative log-likelihood subset gradient:
		lnum = Ab{iblock}' * dothi(:);

		switch arg.precon % denominator
		case 'nuyts'
			tmp = aimax * bel; % Nuyts 98 version for ri=0
			lden = Ab{iblock}' * col(tmp);
		case 'classic'
			% generalization to include nonzero ri:
			tmp = (1 - div0(yi(:,ia).*ri(:,ia), yb.^2)) .* bel;
			lden = Ab{iblock}' * col(tmp .* arg.ai(:,ia));
		case 'fast'
			lden = ldenom * (length(ia) / na);
		end

		x = x - arg.step * (lnum ./ lden); % the update!

		x = max(x, arg.pixmin); % lower bound
		x = min(x, arg.pixmax); % upper bound

		if ~isempty(arg.userfun)
			feval(arg.userfun);
		end
	end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end
end
