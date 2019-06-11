 function xs = epl_os_emdp(x, Gb, yi, ci, ri, R, niter, pixmax, chat)
%function xs = epl_os_emdp(x, Gb, yi, ci, ri, R, niter, pixmax, chat)
% Ordered-subsets version of Alvaro De Pierro's modified EM algorithm
% for penalized log-likelihood image reconstruction from Poisson data.
% Based on Mar. 1995 IEEE T-MI paper: A modified expectation maximization...
%
% see eql_sps_os.m for arguments
%
% Copyright 2002-3-14, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

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
if ~isvar('niter')	| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	| isempty(pixmax),	pixmax = inf;	end
if ~isvar('chat')	| isempty(chat),	chat = false;	end

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

if ~isvar('aj') || isempty(aj)	% \sumi \aij
	aj = Gb' * ci(:) / nblock;
	minmax(aj)
end

% for unpenalized case
if isempty(R)
	pgrad = 0;
	Rdenom = zeros(size(x));
end

%
% loop over iterations
%
xs = zeros(numel(x), niter);
if any(x < 0), error 'need x >= 0', end
x = min(x,pixmax);
xs(:,1) = x;
for iter = 2:niter
	if chat, printf('E-QL-OS-EMDP iteration %d', iter-1), end

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

		x = eql_root(Rdenom, (aj + pgrad - Rdenom .* x) / 2, ...
			x .* eterm);
		x = min(x, pixmax);
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,iter) = x;
end
