  function xs = eml_em(x, G, yi, ci, ri, varargin)
%|function xs = eml_em(x, G, yi, ci, ri, [options])
%|
%| E-ML-EM algorithm for image reconstruction from Poisson emission data
%| see em_fbp.m for model, G, yi, ci, ri
%|
%| in
%|	x	[np 1]		initial image guess (column - see examples)
%|
%| option
%|	'niter'	int		# iterations (default: 1+1)
%|	'isave'	[]		list of iterations to archive
%|				(default: 'last')
%|	'Asum' 	[np 1]		A'1 = G' * c	(optional) column sums
%|	userfun	@		user-defined function handle
%|
%| out
%|	xs	[np nsave]	estimates each (saved) iteration
%|
%| Copyright 1998, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('eml_osem_example')
return
end
if nargin < 3, help(mfilename), error(mfilename), end

if nargin > 5 && isnumeric(varargin{1}) % support old style for now
	warn 'obsolete use of eml_em'
	xs = eml_em_old(x, G, yi, ci, ri, varargin{:});
return
end

arg.niter = 1;
arg.isave = 'last';
arg.Asum = [];
arg.userfun = [];

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);

if ~isvar('ci') || isempty(ci)
	ci = 1;	% ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = 0; % zeros(size(yi(:)));
end
if isempty(arg.Asum)
	if length(ci) == 1
		arg.Asum = sum(G)' * ci;
	else
		arg.Asum = G' * ci;
	end
end

if any(arg.Asum <= 0), error 'Asum must be positive. Adjust mask.', end

eml_check(yi, ci, ri);

if any(x <= 0), error 'need x > 0', end

xs = zeros(numel(x), length(arg.isave));
if any(arg.isave == 0)
	xs(:,find(arg.isave == 0)) = x;
end


%
% loop over iterations
%
for iter = 1:arg.niter
	ticker(mfilename, iter, arg.niter)

	%
	% update:
	% x = x .* (A' * (y ./ (A * x + r))) ./ sum(A)'
	% where A = D(ci) G
	%

	yp = ci .* (G * x) + ri; % predicted measurements
	eterm = G' * (ci .* (yi ./ yp));
	x = x .* eterm ./ arg.Asum;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
	end

	if ~isempty(arg.userfun)
		feval(arg.userfun);
	end
end
