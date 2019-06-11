  function xs = epl_os_emdp(x, Gb, yi, ci, ri, R, varargin)
%|function xs = epl_os_emdp(x, Gb, yi, ci, ri, R, [option])
%| Ordered-subsets version of Alvaro De Pierro's modified EM algorithm
%| for penalized log-likelihood image reconstruction from Poisson data.
%| Based on Mar. 1995 IEEE T-MI paper: A modified expectation maximization...
%|
%| see eql_sps_os.m for arguments
%|
%| option
%|	'pixmax'		maximum value (default inf)
%|	'aj'	[np 1]		sum(A)
%|	'niter'	int		# iterations (default: 1+1)
%|	'isave'	[]		list of iterations to archive
%|				(default: 'last')
%|	'userfun' @		user-defined function handle
%|	'chat'	int		default 0
%|
%| out
%|	xs	[np nsave]	estimates each (saved) iteration
%|
%| Copyright 2002-3-14, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('eml_osem_example')
return
end

if nargin < 3, help(mfilename), error(mfilename), end

if nargin > 5 && isnumeric(varargin{1}) % support old style for now
	warn 'obsolete use of epl_os_emdp'
	xs = epl_os_emdp_old(x, Gb, yi, ci, ri, varargin{:});
return
end

arg.niter = 1;
arg.isave = 'last';
arg.pixmax = inf;
arg.chat = false;
arg.aj = [];

arg.userfun = [];

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

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

if isempty(arg.aj) % \sumi \aij
	arg.aj = Gb' * ci(:) / nblock;
	if arg.chat
		pr minmax(arg.aj)
	end
end

% for unpenalized case
if isempty(R)
	pgrad = 0;
	Rdenom = zeros(size(x));
end

%
% loop over iterations
%
xs = zeros(numel(x), length(arg.isave));
if any(x < 0), error 'need x >= 0', end
x = min(x, arg.pixmax);
if any(arg.isave == 0)
	xs(:,find(arg.isave == 0)) = x;
end

for iter = 1:arg.niter
	if arg.chat, printf('E-QL-OS-EMDP iteration %d', iter-1), end
	ticker(mfilename, iter, arg.niter)

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);

%		ia = iblock:nblock:na;
%		li = Gb{iblock} * x;
%		li = reshape(li, nb, length(ia));
%		yb = ci(:,ia) .* li + ri(:,ia);		% predicted measurements
%		if any(yi(:,ia) & ~yb), error 'model mismatch', end
%		yb(yb == 0) = inf;
%		eterm = Gb{iblock}' * col(ci(:,ia) .* (yi(:,ia) ./ yb));
		eterm = eml_eterm(x, Gb, yi, ci, ri, iblock);

		if ~isempty(R)
			pgrad = R.cgrad(R, x);
			Rdenom = R.denom(R, x);
		end

		x = eql_root(Rdenom, (arg.aj + pgrad - Rdenom .* x) / 2, ...
			x .* eterm);
		x = min(x, arg.pixmax);
	end

	if arg.chat, printf('Range %g %g', min(x), max(x)), end

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end

	if ~isempty(arg.userfun)
		feval(arg.userfun);
	end
end
