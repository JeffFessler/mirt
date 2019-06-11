 function R = Rmem(kappa, varargin)
%function R = Rmem(kappa, [options])
%
% Build roughness penalty regularization "object" based on C = Cdiff() object,
% for regularized solutions to inverse problems.
% This version is somewhat less flexible than Robject() but uses less memory.
%
% General form of nonquadratic penalty function:
% R(x) = \sum_{m=1}^M \beta_m R_m(x),
% where R_m(x) = \sum_k w_k \pot_m([C_m x]_k),
% where [C_m x]_k = \sum_j c_{l,k,j} x_j.
% The matrix C_m is essentially np times np.
% Note that w_k = kappa_k^2 is independent of "m" (offset direction).
%
% Penalty gradient is \sum_m C_m' D_m C_m x,
% where D_m = diag{w_k \wpot_m([C_m x]_k)} and \wpot_m(t) = \pot_m(t) / t.
%
% in
%	kappa	[nx,ny[,nz]]	kappa array, or logical support mask
%
% options
%	'order', 1|2		1st-order or 2nd-order differences (see Cdiff)
%	'offsets', [M] | char
%				offsets to neighboring pixels
%				(see Cdiff for the defaults)
%				use '3d:26' to penalize all 13 pairs of nbrs
%				use '0' for C = I (identity matrix)
%	'beta', [1] | [M]	global regularization parameter(s)
%				default: 2^0
%	'pot_arg', {} 		arguments to potential_func()
%				e.g., {'huber', delta}, or cell{M} array
%				default: {'quad'} for quadratic regularization.
%	'type_denom', ''	type of "denominator"
%				(for quadratic surrogates like SPS)
%		'matlab'	denominator for SPS
%		'aspire'	denominator for SPS that matches aspire
%		'none'		no denominator (default)
%				(since R.E and R.denom only needed for SPS)
%	'distance_power', [1]	See Rweights.m
%	'mask'			Override default: mask = (kappa ~= 0)
%
% out
%	R	strum object with methods: 
%	R.penal(x)	evaluates R(x)
%	R.cgrad(x)	evaluates \cgrad R(x) (column gradient)
%	R.denom(x)	evaluates denominator for separable surrogate
%	[pderiv pcurv] = feval(R.dercurv, R, C1*x) derivatives and curvatures
%			for non-separable parabola surrogates
%	R.diag		diagonal of Hessian of R (at x=0), for preconditioners.
%	R.C1		differencing matrix, with entries 1 and -1,
%			almost always should be used in conjunction with R.wt
%
% Typical use:	mask = true([128 128]); % or something more conformal
%		R = Rmem(mask, 'beta', 2^7);
%
% Copyright 2006-5-23, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(kappa, 'test'), run_mfile_local 'Rmem_test', return, end

% option defaults
R.pot_arg = {'quad'};
R.beta = 2^0;
R.type_denom = 'none';
R.distance_power = 1;
R.order = 1; % 1st-order differences, used in call to Cdiff
R.mask = [];
R.offsets = [];

% parse name/value option pairs
R = vararg_pair(R, varargin);

R.dim = size(kappa);
R.kappa2 = single(kappa(:) .^ 2);

% dimensions, and default offsets
R.offsets = penalty_offsets(R.offsets, size(kappa));
R.offsets = int32(R.offsets);
R.M = length(R.offsets);

% beta setup, including distance_power effect
R.beta = R.beta(:) ./ penalty_distance(R.offsets(:), R.dim) .^ R.distance_power;

% potential setup
if iscell(R.pot_arg{1})
	if length(R.pot_arg) ~= R.M
		error 'pot_arg size mismatch with offsets'
	end
else
	R.pot_arg = { R.pot_arg };
end

R.isquad = true;
for mm=1:R.M
	R.isquad = R.isquad & streq(R.pot_arg{min(mm,end)}, 'quad');
end

% mask
if isempty(R.mask)
	R.mask = kappa ~= 0; % default is to infer from kappas
else
	if ~islogical(R.mask), error 'mask must be logical', end
	R.kappa2 = R.kappa2 .* R.mask(:); % 'mask' the kappas
end
for id=1:ndims(R.mask)-1
	tmp = shiftdim(R.mask, id-1);
	tmp = reshape(tmp, size(tmp,1), []);
	if any(tmp(1,:)) || any(tmp(end,:))
		printm('WARN: mask is nonzero at boundaries along dim=%d', id)
	end
end

%
% build plain sparse differencing object "C1", containing 1 and -1 values
% (or identity matrix)
%
R.C_is_I = any(R.offsets == 0);
if R.C_is_I
	if any(R.offsets), error 'identity is offsets=0', end
	R.C1 = diag_sp(ones(sum(R.mask(:)),1)); % identity matrix

else
	R.C1 = Cdiff(R.mask, 'offsets', R.offsets, ...
			'edge_type', 'none', 'order', R.order);
end
R.np = sum(R.mask(:));

% trick: for quadratic case, provide a R.C object where C = sqrt(W) * C1
R.C = Rmem_setup_C(R);


%
% desired potential functions
%
for mm=1:R.M
	R.pot{mm} = potential_func(R.pot_arg{min(mm,end)}{:});
end

%
% precompute denominator for separable quadratic case if needed
%
if ~streq(R.type_denom, 'none') & R.isquad
	R.denom_max0 = Rmem_denom_max0(R);
end

R.diff_str_forw = sprintf('diff%d,forw%d', R.order, 1);
R.diff_str_back = sprintf('diff%d,back%d', R.order, 1);

% strum methods
% trick: for backwards compatibility, all these *require* that R
% is passed (as dummy argument) even though "strum" does that.
R.dercurv = @Rmem_dercurv; % trick: requires feval()
meth = {...
	'penal', @Rmem_penal, '(R, x)'; ...
	'cgrad', @Rmem_cgrad, '(R, x)'; ...
	'diag', @Rmem_diag, '(R)'; ...
	'denom', @Rmem_denom, '(R, x)'; ...
	'numer_pl_pcg_qs_ls', @Rmem_numer_pl_pcg_qs_ls, '(R, x1, x2)'; ...
	'denom_pl_pcg_qs_ls', @Rmem_denom_pl_pcg_qs_ls, '(R, x1, x2)'; ...
	'numer_denom_pl_pcg_qs_ls', @Rmem_numer_denom_pl_pcg_qs_ls, '(R, x1, x2)'
	};
R = strum(R, meth);


%
% Rmem_penal()
% trick: allow multiple realizations of x
%
function penal = Rmem_penal(R, dummy, x)
if size(x,1) == R.np
	x = embed(x, R.mask);
end
x = reshape(x, prod(R.dim), []);
LL = ncol(x);
penal = zeros(LL,1);
for ll=1:LL
	penal(ll) = Rmem_penal1(R, x(:,ll));
end

%
% Rmem_penal1()
% penalty value for a single image x
%
function penal = Rmem_penal1(R, x)
penal = 0;
for mm=1:R.M
	d = penalty_mex_call(R.diff_str_forw, x, ...
		R.offsets(mm), 2); % 2 because x is a single column
	pot = R.pot{mm};
	d = pot.potk(pot, d);
	penal = penal + R.beta(mm) * sum(R.kappa2 .* d);
end


%
% Rmem_cgrad()
%
function cgrad = Rmem_cgrad(R, dummy, x)
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
	cgrad(:,ll) = Rmem_cgrad1(R, x);
end
if flag_column
	cgrad = cgrad(R.mask(:), :);
end
cgrad = reshape(cgrad, siz);

%
% Rmem_cgrad1()
%
function cgrad = Rmem_cgrad1(R, x)
cgrad = 0;
for mm=1:R.M
	d = penalty_mex_call(R.diff_str_forw, x, ...
		R.offsets(mm), 2); % 2 because a single column
	pot = R.pot{mm};
	d = R.kappa2 .* pot.dpot(pot, d);
	d = penalty_mex_call(R.diff_str_back, d, ...
		R.offsets(mm), 2); % 2 because a single column
	cgrad = cgrad + R.beta(mm) * d;
end
cgrad = cgrad .* R.mask(:);


%
% Rmem_diag()
% evaluate diagonal of Hessian of (quadratic surrogate for) R at x=0
% \sum_k w_k |c_kj|^2 \wpot(0)
%
function rjj = Rmem_diag(R, dummy)
error 'not done, ask jeff'
if R.C_is_I
	error 'not done, ask jeff'
else
	t = reshape(t, [R.dim R.M]);
	rjj = penalty_mex('diff1,back2', t, R.offsets);
	rjj = double(rjj(R.mask));

	for mm=1:R.M
		pot = R.pot{mm};
		t = single(R.beta(mm) * R.kappa2 .* pot{mm}.wpot(pot, 0));
		d = penalty_mex_call(R.diff_str_forw, x, ...
			R.offsets(mm), 2); % 2 because a single column
		d = R.kappa2 .* pot.wpot(pot, d) .* pot.dpot(pot, d);
		d = penalty_mex_call(R.diff_str_back, d, ...
			R.offsets(mm), 2); % 2 because a single column
		cgrad = cgrad + R.beta(mm) * d;
	end
end


%
% Rmem_denom()
% jth penalty separable surrogate curvature is d_j = \sumk |\ckj| \ck \wpotk
% where \ck = \sumj |\ckj|.
% Here, ck = 2 since there is a +1 and a -1 per row of C1 (for 1st-order) 
% Also, |\ckj| = |\ckj|^2 since \ckj = +/- 1, so we can use 'diff1,back2'
%
function denom = Rmem_denom(R, x)
error 'not done, ask jeff'
if streq(R.type_denom, 'none')
	error 'denom not initialized'
end

if R.isquad
	denom = R.denom_max0;
else
	t = single(R.wt .* R.pot.wpot(R.pot, R.C1 * x));
	if R.C_is_I
		denom = t;
	else
		t = reshape(t, [R.dim R.M]);
		if R.order == 1
			t = 2 * t; % because (1,-1)
			denom = penalty_mex('diff1,back2', t, R.offsets);
		elseif R.order == 2 % (-1,2,1) so ck = 4
			t = 4 * t; % 4 = |-1| + |2| + |-1|
			R.denom_max0 = penalty_mex('diff2,backA', ...
					single(t), R.offsets);
			warning 'todo: order=2 not tested'
		else
			error 'order not done'
		end
		denom = denom(R.mask);
	end
end


%
% Rmem_denom_max0()
%
function denom_max0 = Rmem_denom_max0(R)
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
		denom_max0 = R.denom_max0(R.mask);
	end
otherwise
	fail('Unknown type_denom: "%s"', R.type_denom)
end


%
% Rmem_numer_pl_pcg_qs_ls()
%
function numer = Rmem_numer_pl_pcg_qs_ls(R, x1, x2)
numer = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d1 = R.beta(mm) * R.kappa2 .* pot.dpot(pot, d1);
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	numer = numer + sum(col(d1 .* d2));
end


%
% Rmem_denom_pl_pcg_qs_ls()
%
function denom = Rmem_denom_pl_pcg_qs_ls(R, x1, x2)
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d1 = R.beta(mm) * R.kappa2 .* pot.wpot(pot, d1);
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	denom = denom + sum(col(d1 .* (d2.^2)));
end


%
% Rmem_numer_denom_pl_pcg_qs_ls()
%
function out = Rmem_numer_denom_pl_pcg_qs_ls(R, x1, x2)
numer = 0;
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	d1_numer = R.beta(mm) * R.kappa2 .* pot.dpot(pot, d1);
	d1_denom = R.beta(mm) * R.kappa2 .* pot.wpot(pot, d1);
	numer = numer + sum(col(d1_numer .* d2));
	denom = denom + sum(col(d1_denom .* (d2.^2)));
end
out = [numer denom];

%
% Rmem_dercurv()
% evaluate \dpoti and \wpoti
%
function [deriv, curv] = Rmem_dercurv(R, C1x)
error 'not done, ask jeff'
deriv = R.wt .* R.pot.wpot(R.pot, C1x) .* C1x; % fix: use dpot?
curv = R.wt .* R.pot.wpot(R.pot, C1x);



%
% Compute both cgrad and denom (of separable surrogate) efficiently.
% Unused for now, but could be used if inlines are too inefficient
% since both R.cgrad and R.denom use C*x so there is redundancy.
% What we really need is an inline that has two output arguments.
% No, the feval with a function_handle will suffice!
%
%function [cgrad, denom] = Rmem_cgrad_denom(R, x)
%if R.isquad
%	cgrad = R.C' * (R.C * x);
%	denom = R.denom;
%else
%	Cx = R.C * x;
%	wx = R.wt .* R.pot.wpot(R.pot, Cx);
%	cgrad = R.C' * (wx .* Cx);
%	denom = R.E * wx;
%end



%
% Rmem_setup_C()
%
function C = Rmem_setup_C(R)
if R.isquad
%	R.C = R.C1;
%	R.C.cascade_after = diag_sp(sqrt(R.wt));
	R.Cquad_after = Fatrix(size(R.C1,1)*[1 1], R, ...
		'forw', @Rmem_Cquad_after_forw)
	R.C = R.Cquad_after * R.C1;
else
	C = [];
end

%
% Rmem_Cquad_after_forw()
% kron(sqrt(beta_m), diag(kappa))
%
function d = Rmem_Cquad_after_forw(R, d)
tmp = R.kappa2 * R.beta(:)'; % [np,M]
d = repmat(sqrt(tmp(:)), 1, ncol(d)) .* d;
