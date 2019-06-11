 function inv2 = de_ftab_inv2(fit, sl, varargin)
%function inv2 = de_ftab_inv2(fit, sl, [options])
%|
%| Build object that does polynomial inverse of BH function.
%| (To map log data into corresponding to line-integrals of material density.)
%|
%| in
%|	fit	strum	initialized by de_ftab_fit()
%|	sl	cell	{s1, ..., sL} see sls.sl from de_ftab_sls.m
%|
%| option
%|	'itype'	char	'T' or 'fixlinear' or 'fit' todo
%|				'fixlinear' is default which means
%|				use ftab.T for linear terms and then
%|				fit higher-order polynomial for rest
%|	'T'	[L M]	ftab.T - linearization at 0, if desired
%|	'show'	0|1	plot it?
%|
%| out
%|	inv2	strum
%|	methods:
%|		inv2.fun(hf) 		map log values into corrected values
%|		inv2.plot(fit, sl)	show fits
%|
%| Copyright 2008-09-28, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(fit, 'test'), fail('run de_ftab test'), return, end

% option defaults
arg.emax = 4;
arg.T = [];
arg.esum = 6;
%arg.itype = 'T';
%arg.itype = 'fit'; % general fit
arg.itype = 'fixlinear'; % fix linear coefficients based on T, fit rest
%arg.wt = sll{2}.^3;	% weighting
%arg.wt = 1 ./ (0*sll{2} + 1);	% weighting
arg.show = false;

arg = vararg_pair(arg, varargin);

inv2 = de_ftab_inv2_setup(fit, sl, arg);


%
% de_ftab_inv2_setup()
% add methods to ftab related to polynomial approximation to its inverse
%
function inv2 = de_ftab_inv2_setup(fit, sl, arg)

LL = fit.LL;
MM = fit.MM;

if fit.MM < fit.LL
	warn('not making inv2() for MM=%d < LL=%d', fit.MM, fit.LL)
	inv2 = [];
return
end

switch arg.itype
case 'T' % 1st-order poly based on T
	if isempty(arg.T), fail 'need T', end
	coef = arg.T'; % [K L], K=M
	expo = eye(MM); % [K M], K=M

case {'fit', 'fixlinear'}
	sll = ndgrid_jf('mat', sl{:}); % [(Ns) L]
	sll = reshapee(sll, [], LL); % [*Ns L]
	ftmp = fit.fmfun(sll); % [*Ns M]

	expo = de_poly_expo(MM, 'emax', arg.emax, 'esum', arg.esum);
	expo(1,:) = []; % remove DC term
	basis = de_poly_eval(ftmp, [], expo, 'basis', 1); % [*Ns K]

% todo: consider weighting?
	if streq(arg.itype, 'fixlinear')
		if isempty(arg.T), fail 'need T', end
		jf_equal(expo(1:MM,1:MM), eye(MM)) % linear terms
		KK = size(expo,1); % # of coefficients
		coef = zeros(KK, LL);
		coef(1:MM,:) = arg.T'; % [M L]
		lin = de_poly_eval(ftmp, coef(1:MM,:), expo(1:MM,:)); % [*Ns L]
		basis = basis(:,(MM+1):end); % [*Ns K-M]
		coef((MM+1):end,:) = basis \ (sll - lin); % [K-M L]
	else
		coef = basis \ sll; % [K L]
	end

otherwise
	fail('bad itype %s', arg.itype)
end

arg.coef = coef;
arg.expo = expo;

if 1 % display polynomial fit errors
	tmp = de_poly_eval(ftmp, coef, expo);
	for ll=1:LL
		printm('worst inverse poly error l=%d: %g of %g', ll, ...
			max(abs(tmp(:,ll) - sll(:,ll))), max(sll(:,ll)))
	end
end

%
% old way: (todo: remove?)
% Build polynomial approximation to shat = (T F)^{-1}(z).
% Use weighted fitting so that soft tissue is fit well!
%

if 0	% this polynomial way seems insufficiently accurate!?
	ftab.inv.basis_func = ir_poly2_fun(3);
	tmp = ftab.inv.basis_func(z1(:), z2(:));
	ftab.inv.nbasis = length(ftab.inv.basis_func(0,0));
	ftab.inv.coef = zeros(ftab.inv.nbasis,2);
	for ll=1:ftab.LL
		ftab.inv.coef(:,ll) = (diag_sp(wt(:)) * tmp) \ ...
			(wt(:) .* sll{ll}(:));
	end
end

%ftab.inv.eval = @(ftab,fha) de_ftab_invert(ftab,fhat);

fun = @(inv2, fhat) de_ftab_inv2_eval(inv2, fhat);
meth = {'fun', fun, '(fhat [() M]) -> [() L]';
	'plot', @de_ftab_inv2_plot, '(fit, sl)'};
inv2 = strum(arg, meth);



%
% de_ftab_inv2_eval()
% todo
%
function shat = de_ftab_inv2_eval(inv2, fhat)
shat = de_poly_eval(fhat, inv2.coef, inv2.expo);


%
% de_ftab_inv2_plot()
% see how well a grid of s values is recovered
%
function out = de_ftab_inv2_plot(inv2, fit, sl)
LL = fit.LL;
if LL == 2
	sll = ndgrid_jf('mat', sl{:});
	ftmp = fit.fmfun(sll);
	shat = inv2.fun(ftmp);
	if im
		clf
	%	plot(sl{1}, sl{2}, shat, '.')
		plot(shat(:,:,1), shat(:,:,2), '.')
		axis tight
		xlabel 's1', ylabel 's2'
		title 's -> (fit) -> f -> (inv2) -> s'
	prompt
	end
else
	warn 'inv2_plot only for L=2'
end

if nargout, out = []; end
