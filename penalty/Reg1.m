 function R = Reg1(kappa, varargin)
%function R = Reg1(kappa, [options])
%|
%| Build roughness penalty regularization "object" based on Cdiff1() objects,
%| for regularized solutions to inverse problems.
%| This version supercedes Robject() by providing its capabilities while
%| also providing options that use less memory.  By default it tries to use
%| a mex version (penalty_mex.mex*) but if that fails it reverts to a pure
%| matlab form that should be completely portable.
%|
%| General form of (possibly nonquadratic) penalty function:
%|	R(x) = \sum_{m=1}^M sum_n w[n;m] potential_m( [C_m x]_n )
%| where M is the number of neighbors (offsets), and n=1,...,N=numel(x)
%| and C_m is a (N x N) differencing matrix in the direction off_m.
%|
%| See Rweights() for w[n;m] options, only some of which use low memory.
%|
%| Penalty gradient is \sum_m C_m' D_m(x) C_m x,
%| where D_m(x) = diag{w[n;m] \wpot_[n;m]([C_m x]_n)} and \wpot(t) = \pot(t) / t
%|
%| in
%|	kappa	[(N)]		kappa array, or logical support mask
%|
%| options
%|	'type_penal'		'def' | '' : try 'mex', otherwise 'mat'
%|				'mex' : compute penalty gradient with mex call
%|				'mat' : compute penalty gradient via Cdiff1
%|				'zxy' : mex penalty for zxy ordered 3D image
%|					(caution: see zxy conventions below)
%|	'type_diff'		'def|ind|mex|sparse|...' (see Cdiff1)
%|	'order', 1|2		1st-order or 2nd-order differences (see Cdiff1)
%|	'offsets', [M] | char
%|				offsets to neighboring pixels
%|					(see Cdiff1 for the defaults)
%|				use '3d:26' to penalize all 13 pairs of nbrs
%|				use "0" for C = I (identity matrix)
%|	'beta', [1] | [M]	global regularization parameter(s)
%|				default: 2^0
%|	'pot_arg', {} 		arguments to potential_fun()
%|					e.g., {'huber', delta}, or cell{M} array
%|				default: {'quad'} for quadratic regularization.
%| ?	'pre_denom_sqs1_x0'	precompute denominator for SQS at x=0? (def: 0)
%|	'type_denom', ''	type of "denominator"
%|					(for quadratic surrogates like SPS)
%|					todo: improve documentation!
%|		'matlab'	denominator for SPS
%|					todo: precompute?  or just on the fly?
%|		'aspire'	denominator for SPS that matches aspire
%|		'none'		no denominator precomputation (default)
%|	'distance_power', 0|1|2	See Rweights.m
%|	'user_wt', [(N) M]	""
%| ?	wt_use_mex 0|1		""
%|	'nthread'	(int)	# of threads (for type_penal='mex' only)
%|	'mask'			Override default: mask = (kappa ~= 0)
%|				(Only use if you sure know what you're doing!)
%|
%| out
%|	R	strum object with methods: 
%|	R.penal(x)	evaluates R(x)
%|	R.cgrad(x)	evaluates \cgrad R(x) (column gradient)
%|	R.denom_sqs1(x)	evaluates denominator for separable quadratic surrogate
%|				\sum_m |C_m|' |D_m(x)| (|C_m| 1)
%|	R.denom(x)	evaluates denominator for separable surrogate
%|	[pderiv pcurv] = feval(R.dercurv, R, R.C1*x) derivatives and curvatures
%|				for non-separable parabola surrogates
%|	R.diag		diagonal of Hessian of R (at x=0), for preconditioners.
%|	R.C1		finite differencing matrix; for 1st-order differences
%|			nonzero entries of each row are 1 and -1;
%|			for 2nd-order differences each row contains {-1 2 -1}.
%|			Almost always should be used in conjunction with R.wt.
%|	R.C		diag(sqrt(w)) * C1, only for quadratic case!
%|	R.wt		see Rweights
%|
%| Typical use:	mask = true(128, 128); % or something more conformal
%|		R = Reg1(mask, 'beta', 2^7);
%|
%| zxy conventions:
%|	kappa must be given as [nz nx ny] by user
%|	offsets must be given w.r.t. zxy; use reg_offset_xyz_to_zxy() if needed.
%|	To ensure the user is aware of this, must set flag 'offsets_is_zxy'
%|
%| Note: R.denom() is slightly larger near image borders (but still valid) 
%| when using the mex version compared to the mat version in the formula above.
%|
%| Copyright 2006-12-6, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(kappa, 'test'), run_mfile_local 'Reg1_test', return, end

% option defaults
arg.type_penal = '';
%arg.edge_type = 'simple'; % saves memory
arg.edge_type = 'tight'; % because mex uses this
arg.type_diff = ''; % defer to Cdiff1 default
arg.pot_arg = {'quad'};
arg.beta = 2^0;
arg.type_denom = 'none';
arg.pre_denom_sqs1_x0 = false;
arg.distance_power = 1; % perhaps '2' would be better
arg.wt_use_mex = false;
arg.user_wt = []; % for Rweights
arg.order = 1; % 1st-order differences
arg.mask = [];
arg.offsets = [];
arg.offsets_is_zxy = false;
arg.control = []; % todo: document this
arg.type_wt = '';
arg.nthread = 1;

if isempty(arg.type_wt)
	if numel(kappa) <= 128^2
		arg.type_wt = 'pre'; % small enough to precompute
	else
		arg.type_wt = 'fly'; % saves memory
	end
end

% parse name/value option pairs
arg = vararg_pair(arg, varargin);

% default offsets
if ischar(arg.offsets) && streq(arg.type_penal, 'zxy')
	arg.offsets_is_zxy = true; % trick
end
arg.offsets = penalty_offsets(arg.offsets, size(kappa));
arg.M = numel(arg.offsets);

if prod(size(kappa)) < arg.order + 1
	warn('image size %d too small for order=%d', prod(size(kappa)), arg.order)
%{
	warn('image size %d too small for order=%d, forming null penalty', ...
		prod(size(kappa)), arg.order)
	R = Reg1_null(size(kappa));
%	arg.order = max(prod(size(kappa))-1, 0);
	arg.order = 0;
	warn('reverting to order=%d', arg.order)
	arg.offsets = 0;
	arg.M = 1;
%}
end

if isempty(arg.control)
	switch arg.order
	case 0
		arg.control = 0; % irrelevant
	case 1
		arg.control = 1; % irrelevant
	case 2
		if isempty(arg.user_wt)
			arg.control = 2; % slightly faster
		else
			arg.control = 1; % control=2 not done for user_wt
		end
	otherwise
		fail('order=%d control default not done', arg.order)
	end
end

% dimensions
arg.dim = size(kappa);
if ndims(kappa) == 2 && size(kappa,2) == 1 % 1D 2015-08-17
	arg.dim = size(kappa, 1);
end

R.kappa = kappa; % save kappa for later reference if needed

if ~isreal(kappa), fail 'kappa must be real', end
if any(kappa(:) < 0), warn 'kappa has negative values?', end

% potential function setup
if iscell(arg.pot_arg{1})
	if length(arg.pot_arg) ~= arg.M
		error 'pot_arg size mismatch with offsets'
	end
else
	arg.pot_arg = { arg.pot_arg };
end

% determine if quadratic
arg.isquad = true;
for mm=1:arg.M
	arg.isquad = arg.isquad & streq(arg.pot_arg{min(mm,end)}{1}, 'quad');
end

% mask
if isempty(arg.mask)
	arg.mask = kappa ~= 0; % default is to infer from kappas
else
	if ~islogical(arg.mask), error 'mask must be logical', end
end
arg.np = sum(arg.mask(:));
%mask_border_check(arg.mask);

% weights
arg.wt = Rweights(kappa, arg.offsets, 'type_wt', arg.type_wt, ...
		'edge_type', arg.edge_type, 'beta', arg.beta, ...
		'order', arg.order, 'distance_power', arg.distance_power, ...
		'user_wt', arg.user_wt, ...
		'use_mex', 0); % todo: because 'simple'

% differencing objects
arg.C1s = cell(arg.M,1);
for mm=1:arg.M
	arg.C1s{mm} = Cdiff1(arg.dim, 'type_diff', arg.type_diff, ...
		'offset', arg.offsets(mm), 'order', arg.order);
end

% desired potential function handles
for mm=1:arg.M
	tmp = arg.pot_arg{min(mm,end)};
	arg.pot{mm} = potential_fun(tmp{:});
end

% here the different types of penalty implementations depart

%if isempty(arg.type_penal) && ~isempty(arg.user_wt)
%	warn '"mex" does not support "user_wt" so using "mat"'
%	arg.type_penal = 'mat';
%end

if isempty(arg.type_penal) || streq(arg.type_penal, 'def')
	if exist('penalty_mex') == 3
		arg.type_penal = 'mex'; % default: try mex if available
	else
		arg.type_penal = 'mat';
	end
end

if xor(streq(arg.type_penal, 'zxy'), arg.offsets_is_zxy)
	fail 'offsets and kappa must be zxy iff type_penal is zxy'
end

switch arg.type_penal
case 'mat'
	if arg.nthread ~= 1, warn 'nthread != 1 ignored', end
	R = Reg1_setup_mat(arg, kappa);

case 'mex'
	% the 'mex' code has various restrictions (for efficiency)
	% so here we make sure that those limitations are obeyed.
	if ~streq(arg.type_wt, 'pre') && ~streq(arg.type_wt, 'fly')
		warn 'for "mex" type_penal, "type_wt" should be "pre" or "fly"'
		warn('are you sure about what you are doing? (%s)', arg.type_wt)
	end
% todo	if ~streq(arg.edge_type, 'simple')
	if ~streq(arg.edge_type, 'tight')
		fail 'for "mex" penalty type, "edge_type" should be "tight"'
	end
	if streq(arg.type_diff, 'sparse') || streq(arg.type_diff, 'spmat')
		warn 'for "mex" penalty type, "type_diff" should be "ind|mex"'
		warning 'are you sure you know what you are doing?'
	end
	if length(arg.pot_arg) > 1
		fail 'for "mex" penalty type, pot_arg must be only cell(1)'
	end

	R = Reg1_setup_mex(arg, kappa);

case 'zxy'
	if ~isempty(arg.user_wt)
		fail '"zxy" type_penal does not support "user_wt" option'
	end
	R = Reg1_setup_zxy(arg, kappa);

otherwise
	fail('bad penalty type "%s"', arg.type_penal)
end


% Reg1_setup_mat()
% matlab-based calculations of penalty gradient.
% May be partially mex-based depending on Cdiff1() options, but even then
% the Cdiff1 objects are the only mex aspect; all else is matlab.
function R = Reg1_setup_mat(arg, kappa)

%arg.kappa2 = arg.kappa2 .* arg.mask(:); % 'mask' the kappas
%arg.kappa2 = single(kappa(:) .^ 2);

% todo: identity?
%arg.C_is_I = isequal(arg.offsets, [0]);

% precompute denominator for separable quadratic surrogate if requested
if arg.pre_denom_sqs1_x0
	arg.denom_sqs1_x0 = Reg_mat_denom_max0(arg);
end

% strum methods
% trick: for backwards compatibility, all these *require* that R
% is passed (as dummy argument) even though "strum" does that.
arg.cgrad1_fun = @Reg1_mat_cgrad1; % used in Reg1_com_cgrad
arg.dercurv = @Reg1_com_dercurv; % trick: requires feval(), in private/
meth = {...
	'C1', @Reg1_com_C1, '(R)'; ...
	'C', @Reg1_com_C, '(R)'; ...
	'penal', @Reg1_com_penal, '(R, x)'; ...
	'cgrad', @Reg1_com_cgrad, '(R, x)'; ...
	'egrad', @Reg1_com_egrad, '(R, x, delta)'; ...
	'denom_sqs1', @Reg1_mat_denom_sqs1, '(R, x)'; ...
	'denom', @Reg1_mat_denom, '(R, x)'; ...
%	'diag', @Reg1_mat_diag, '(R)'; ...
%	'numer_pl_pcg_qs_ls', @Reg1_mat_numer_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'denom_pl_pcg_qs_ls', @Reg1_mat_denom_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'numer_denom_pl_pcg_qs_ls', @Reg1_mat_numer_denom_pl_pcg_qs_ls, '(R, x1, x2)'
	};
R = strum(arg, meth);


%
% 'mat' methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


% Reg1_mat_cgrad1()
% this is called by Reg1_com_cgrad()
function cgrad = Reg1_mat_cgrad1(sr, x)
cgrad = 0;
for mm=1:sr.M
	d = sr.C1s{mm} * x;
	pot = sr.pot{mm};
	wt = pot.wpot(d);
	wt = wt .* sr.wt.col(mm);
	tmp = sr.C1s{mm}' * (wt .* d);
	cgrad = cgrad + tmp;
end
cgrad = cgrad .* sr.mask(:);


% Reg1_mat_denom()
function denom = Reg1_mat_denom(sr, dummy, x)
switch sr.type_denom
case {'none', 'matlab'}
	;
case {'aspire'}
	warn 'aspire denom not done, using matlab'
otherwise
	fail('unknown type_denom %s', sr.type_denom)
end
%if sr.isquad
%	denom = R.denom_max0;
%else
denom = Reg1_mat_denom_sqs1(sr, dummy, x);
%end


% Reg1_mat_denom_sqs1()
% penalty "denominator" term (curvature) for separable quadratic surrogate.
% d_j = \sumk |\ckj| \ck \wpot([Cx]_k), where \ck = \sumj |\ckj|.
% use x=0 (default) to get "maximum" denominator
function denom = Reg1_mat_denom_sqs1(sr, dummy, x)

if ~isvar('x') || isempty(x)
	x = zeros(size(sr.mask));
end
denom = Reg1_mat_denom_sqs1_cell(sr.C1s, sr.mask, sr.pot, sr.wt, x);


% Reg1_setup_mex()
% mex-based calculations of penalty gradient.
function R = Reg1_setup_mex(arg, kappa)

% penalty setup
arg.pot_type = arg.pot_arg{1}{1};
arg.pot_params = cat(2, arg.pot_arg{1}{2:end});
dis = penalty_distance(arg.offsets(:), arg.dim) .^ arg.distance_power;
arg.beta = arg.beta(:) ./ dis;

beta_arg = single(arg.beta);

if ~isempty(arg.user_wt) % 2013-01-26 added support for user_wt
	tmp = size(arg.user_wt);
	if tmp(end) ~= arg.M
		fail('user_wt must be [(N) M]')
	end
	tmp = reshapee(arg.user_wt, [], arg.M); % [nn M]
	tmp = repmat(arg.beta, [1 nrow(tmp)]) .* tmp'; % [M nn] is_beta_array!
	beta_arg = single(tmp);
end

% arguments that will be passed to penalty_mex() later
arg.cgrad1_str = 'cgrad,offset';
arg.denom_sqs1_str = 'denom,offset';
arg.penal_str = 'penal,offset';
arg.cdp_arg = { single(kappa), int32(arg.offsets), ...
		beta_arg, arg.pot_type, single(arg.pot_params) ...
		int32(arg.order), int32(arg.control), int32(arg.nthread) };

%arg.diff_str_forw = sprintf('diff%d,forw%d', arg.order, 1);
%arg.diff_str_back = sprintf('diff%d,back%d', arg.order, 1);

% strum methods
% trick: for backwards compatibility, all these *require* that R
% is passed (as dummy argument) even though "strum" does that.
arg.cgrad1_fun = @Reg1_mex_cgrad1; % used in Reg1_com_cgrad
arg.dercurv = @Reg1_com_dercurv; % trick: requires feval(), in private/
meth = {...
	'C1', @Reg1_com_C1, '(R)'; ...
	'C', @Reg1_com_C, '(R)'; ...
	'penal', @Reg1_com_penal, '(R, x)'; ...
	'cgrad', @Reg1_com_cgrad, '(R, x)'; ...
	'egrad', @Reg1_com_egrad, '(R, x, delta)'; ...
	'denom_sqs1', @Reg1_mex_denom_sqs1, '(R, x)'; ...
	'denom', @Reg1_mex_denom, '(R, x)'; ...
%	'diag', @Reg1_diag, '(R)'; ...
%	'numer_pl_pcg_qs_ls', @Reg1_numer_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'denom_pl_pcg_qs_ls', @Reg1_denom_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'numer_denom_pl_pcg_qs_ls', @Reg1_numer_denom_pl_pcg_qs_ls, '(R, x1, x2)'
	};
R = strum(arg, meth);


%
% 'mex' methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


% Reg1_mex_cgrad1()
% this is called by Reg1_com_cgrad()
function cgrad = Reg1_mex_cgrad1(sr, x)
tmp = sr.cdp_arg; % trick: strum limitation
cgrad = penalty_mex(sr.cgrad1_str, single(x), tmp{:});
%cgrad = single(cgrad);


% Reg1_mex_denom()
function denom = Reg1_mex_denom(sr, dummy, x)
switch sr.type_denom
case {'none', 'matlab'}
	;
case {'aspire'}
	warn 'aspire denom not done, using matlab'
otherwise
	fail('unknown type_denom %s', sr.type_denom)
end
%if sr.isquad
%	denom = R.denom_max0;
%else
denom = Reg1_mex_denom_sqs1(sr, dummy, x);
%end


% Reg1_mex_denom_sqs1()
function denom = Reg1_mex_denom_sqs1(sr, dummy, x)
if ~isvar('x') || isempty(x)
	x = zeros(size(sr.mask));
end
[x ei] = embed_in(x, sr.mask);
tmp = sr.cdp_arg; % trick: strum limitation
denom = penalty_mex(sr.denom_sqs1_str, single(x), tmp{:});
%denom = single(denom);
if ei.column
	denom = denom(sr.mask);
end



%
% common methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Reg1_com_penal()
% trick: allow multiple realizations of x
% x can be col or array
function penal = Reg1_com_penal(sr, dummy, x)
if size(x,1) == sr.np
	x = embed(x, sr.mask);
end
x = reshape(x, prod(sr.dim), []); % [*N L]
LL = ncol(x);
penal = zeros(LL,1);
for ll=1:LL
	penal(ll) = Reg1_com_penal1(sr, x(:,ll));
end


% Reg1_com_penal1()
% penalty value for a single image x (an array stretched into column [*N])
function penal = Reg1_com_penal1(sr, x)
if sr.type_penal == 'mex'
	if ~isreal(x)
		warn 'todo: no mex penalty for complex; revert to matlab'
		sr.type_penal = 'mat'; % trick: revert to matlab
		penal = Reg1_com_penal1(sr, x);
	return
	end
	tmp = sr.cdp_arg; % trick: strum limitation
	penal = penalty_mex(sr.penal_str, single(x), tmp{:});
else
	penal = 0;
	for mm=1:sr.M
		d = sr.C1s{mm} * x;
		pot = sr.pot{mm};
		d = pot.potk(d);
		wt = sr.wt.col(mm);
		% trick: double below helps mat and mex versions match better
		penal = penal + sum(wt .* d, 'double');
	end
end


% Reg1_com_C1()
function C1 = Reg1_com_C1(sr)

C1 = Cdiffs(sr.dim, 'type_diff', sr.type_diff, ...
	'offsets', sr.offsets, 'order', sr.order, 'mask', sr.mask);


% Reg1_com_C()
% diag(sqrt(w)) * C1
function C = Reg1_com_C(sr)

if ~sr.isquad, fail 'C method meaningful only for quadratic potential', end
wmn = zeros(prod(sr.dim), sr.M); % [*N M]
for mm=1:sr.M
	wmn(:,mm) = sr.wt.col(mm);
end
switch class(sr.C1)
case 'fatrix2'
	wmn = reshape(wmn, [sr.dim sr.M]);
	Wh = Gdiag(sqrt(wmn));
case 'Fatrix'
	Wh = diag_sp(sqrt(wmn(:)));
otherwise
	if issparse(sr.C1)
		Wh = spdiag(sqrt(wmn(:)), 'nowarn');
	else
		fail('bug: unknown class "%s"', class(sr.C1))
	end
end
C = Wh * sr.C1; % cascade


% Reg1_com_cgrad()
% x and cgrad are *both* either [(N) (L)] or [np (L)]
function cgrad = Reg1_com_cgrad(sr, dummy, x)

siz = size(x);
[x, ei] = embed_in(x, sr.mask, sr.np);
x = reshape(x, prod(sr.dim), []); % [*N *L]
LL = ncol(x);

cgrad = zeros(size(x), 'single'); % [*N *L]
for ll=1:LL
	cgrad(:,ll) = sr.cgrad1_fun(sr, x(:,ll));
end

if ei.column
	cgrad = cgrad(sr.mask(:), :); % [np *L]
	if LL > 1
		cgrad = ei.shape(cgrad);
	end
else
	cgrad = reshape(cgrad, siz);
end


% Reg1_com_egrad()
% empirical gradient (finite differences) for testing _cgrad
function egrad = Reg1_com_egrad(sr, dummy, x, delta)

egrad = zeros(size(x));
ej = zeros(size(x));
ticker reset
for jj=1:numel(x)
	ticker(mfilename, jj, numel(x))
	ej(jj) = 1;
	egrad(jj) = ( sr.penal([], x + delta * ej) - sr.penal([], x) ) / delta;
	ej(jj) = 0;
end


%
% 'old' methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%{

% Reg1_old_penal()
% trick: allow multiple realizations of x
function penal = Reg1_old_penal(R, dummy, x)
if size(x,1) == R.np
	x = embed(x, R.mask);
end
x = reshape(x, prod(R.dim), []);
LL = ncol(x);
penal = zeros(LL,1);
for ll=1:LL
	penal(ll) = Reg1_old_penal1(R, x(:,ll));
end


% Reg1_old_penal1()
% penalty value for a single image x
function penal = Reg1_old_penal1(R, x)
penal = 0;
for mm=1:R.M
	d = penalty_mex_call(R.diff_str_forw, x, ...
		R.offsets(mm), 2); % 2 because x is a single column
	pot = R.pot{mm};
	d = pot.potk(d);
	penal = penal + R.beta(mm) * sum(R.kappa2 .* d, 'double');
end


% Reg1_old_cgrad()
function cgrad = Reg1_old_cgrad(R, dummy, x)
siz = size(x);
flag_column = 0;
if size(x,1) == R.np
	x = embed(x, R.mask);
	flag_column = 1;
end
x = reshape(x, prod(R.dim), []);
LL = ncol(x);
cgrad = zeros(prod(R.dim),LL, 'single');
for ll=1:LL
	cgrad(:,ll) = Reg1_old_cgrad1(R, x);
end
if flag_column
	cgrad = cgrad(R.mask(:), :);
end
cgrad = reshape(cgrad, siz);


% Reg1_old_cgrad1()
function cgrad = Reg1_old_cgrad1(R, x)
cgrad = 0;
for mm=1:R.M
	d = penalty_mex_call(R.diff_str_forw, x, ...
		R.offsets(mm), 2); % 2 because a single column
	pot = R.pot{mm};
	d = R.kappa2 .* pot.dpot(d);
	d = penalty_mex_call(R.diff_str_back, d, ...
		R.offsets(mm), 2); % 2 because a single column
	cgrad = cgrad + R.beta(mm) * d;
end


% Reg1_old_diag()
% evaluate diagonal of Hessian of (quadratic surrogate for) R at x=0
% \sum_k w_k |c_kj|^2 \wpot(0)
function rjj = Reg1_old_diag(R, dummy)
error 'not done, ask jeff'
if R.C_is_I
	error 'not done, ask jeff'
else
	t = reshape(t, [R.dim R.M]);
	rjj = penalty_mex('diff1,back2', t, R.offsets);
	rjj = rjj(R.mask);

	for mm=1:R.M
		pot = R.pot{mm};
		t = single(R.beta(mm) * R.kappa2 .* pot{mm}.wpot(0));
		d = penalty_mex_call(R.diff_str_forw, x, ...
			R.offsets(mm), 2); % 2 because a single column
		d = R.kappa2 .* pot.wpot(d) .* pot.dpot(d);
		d = penalty_mex_call(R.diff_str_back, d, ...
			R.offsets(mm), 2); % 2 because a single column
		cgrad = cgrad + R.beta(mm) * d;
	end
end


% Reg1_old_denom()
% jth penalty separable surrogate curvature is d_j = \sumk |\ckj| \ck \wpotk
% where \ck = \sumj |\ckj|.
% Here, ck = 2 because there is a +1 and a -1 per row of C1 (for 1st-order) 
% Also, |\ckj| = |\ckj|^2 since \ckj = +/- 1, so we can use 'diff1,back2'
function denom = Reg1_old_denom(R, x)
error 'not done, ask jeff'
if streq(R.type_denom, 'none')
	error 'denom not initialized'
end

if R.isquad
	denom = R.denom_max0;
else
	t = single(R.wt .* R.pot.wpot(R.pot, R.C1 * x));

	if R.C_is_I
		denom = single(t);
	else
		t = reshape(t, [R.dim R.M]);
		if R.order == 1
			t = 2 * t; % because (1,-1)
			denom = penalty_mex('diff1,back2', t, R.offsets);
		elseif R.order == 2 % (-1,2,1) so ck = 4
			t = 4 * t; % 4 = |-1| + |2| + |-1|
			denom = penalty_mex('diff2,backA', ...
					single(t), R.offsets);
			warning 'todo: order=2 not tested'
		else
			error 'order not done'
		end
		denom = single(denom(R.mask));
	end
end


% Reg1_denom_max0()
function denom_max0 = Reg1_denom_max0(R)
switch R.type_denom
case {'matlab', 'aspire'}
	error 'not done due to R.wt'
	t = R.wt .* R.pot.wpot(R.pot, 0);
	if R.C_is_I
		denom_max0 = t;
	else
		t = reshape(t, [R.dim R.M]);
		if R.order == 1 % fix: order=2 denom?
			t = 2 * t; % "2" because (1,-1) differences
			denom_max0 = penalty_mex('diff1,back2', ...
					single(t), R.offsets);
		elseif R.order == 2 % (-1,2,1) so ?
			t = 4 * t; % 4 = |-1| + |2| + |-1|
			denom_max0 = penalty_mex('diff2,backA', ...
					single(t), R.offsets);
			warning 'todo: order=2 not tested'
		else
			error 'order not done for denom'
		end
		denom_max0 = single(R.denom_max0(R.mask));
	end
otherwise
	error('Unknown type_denom: "%s"', R.type_denom)
end


% Reg1_numer_pl_pcg_qs_ls()
function numer = Reg1_numer_pl_pcg_qs_ls(R, x1, x2)
numer = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d1 = R.beta(mm) * R.kappa2 .* pot.dpot(d1);
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	numer = numer + sum(col(d1 .* d2), 'double');
end


% Reg1_denom_pl_pcg_qs_ls()
function denom = Reg1_denom_pl_pcg_qs_ls(R, x1, x2)
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d1 = R.beta(mm) * R.kappa2 .* pot.wpot(d1);
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	denom = denom + sum(col(d1 .* (d2.^2)), 'double');
end


% Reg1_numer_denom_pl_pcg_qs_ls()
function out = Reg1_numer_denom_pl_pcg_qs_ls(R, x1, x2)
numer = 0;
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	d1_numer = R.beta(mm) * R.kappa2 .* pot.dpot(d1);
	d1_denom = R.beta(mm) * R.kappa2 .* pot.wpot(d1);
	numer = numer + sum(col(d1_numer .* d2), 'double');
	denom = denom + sum(col(d1_denom .* (d2.^2)), 'double');
end
out = [numer denom];


% Compute both cgrad and denom (of separable surrogate) efficiently.
% Unused for now, but could be used if anonymous functions are too inefficient
% since both R.cgrad and R.denom use C*x so there is redundancy.
% What we really need is an anonymous function that has two output arguments.
% No, the feval with a function_handle will suffice!
%
%function [cgrad, denom] = Reg1_cgrad_denom(R, x)
%if R.isquad
%	cgrad = R.C' * (R.C * x);
%	denom = R.denom;
%else
%	Cx = R.C * x;
%	wx = R.wt .* R.pot.wpot(R.pot, Cx);
%	cgrad = R.C' * (wx .* Cx);
%	denom = R.E * wx;
%end

% Reg1_null()
% null regularizer where R(x) = 0
function R = Reg1_null(idim)

fail 'todo'

arg.dercurv = @Reg1_com_dercurv; % trick: requires feval()
meth = {...
	'C1', @Reg1_com_C1, '(R)'; ...
	'C', @Reg1_com_C, '(R)'; ...
	'penal', @Reg1_com_penal, '(R, x)'; ...
	'cgrad', @Reg1_com_cgrad, '(R, x)'; ...
	'egrad', @Reg1_com_egrad, '(R, x, delta)'; ...
	'denom_sqs1', @Reg1_mat_denom_sqs1, '(R, x)'; ...
	'denom', @Reg1_mat_denom, '(R, x)'; ...
%	'diag', @Reg1_mat_diag, '(R)'; ...
%	'numer_pl_pcg_qs_ls', @Reg1_mat_numer_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'denom_pl_pcg_qs_ls', @Reg1_mat_denom_pl_pcg_qs_ls, '(R, x1, x2)'; ...
%	'numer_denom_pl_pcg_qs_ls', @Reg1_mat_numer_denom_pl_pcg_qs_ls, '(R, x1, x2)'
	};
R = strum(arg, meth);

function [pderiv pcurv] = Reg1_null_dercurv(R, d)
pderiv = 0;
pcurv = 0;

%}
