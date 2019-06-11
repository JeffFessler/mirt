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
%				(its primary use for consistency with ASPIRE)
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
%	'distance_power', [1]	1 classical (default), 2 possibly improved
%				see penalty_mex()
%				if nonzero, then beta's are scaled
%				by 1 ./ \norm{offset_vector_m}^power
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
R.use_integ = false;

% parse name/value option pairs
R = vararg_pair(R, varargin);

R.dim = size(kappa);


% dimensions, and default offsets
R.offsets = Rmem_offsets(R.offsets, kappa);
R.M = length(R.offsets);
R.use_integ = logical(R.use_integ);
R.kappa = kappa;

% beta setup
if isscalar(R.beta)
	R.beta = R.beta * ones(1, R.M);
elseif ~isequal(size(R.beta), size(R.offsets))
	error 'beta size mismatch with offsets'
end

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
R.kappa2 = single(kappa(:) .^ 2);
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
if ~R.use_integ
  meth = {...
      'penal', @Rmem_penal, '(R, x)'; ...
      'cgrad', @Rmem_cgrad, '(R, x)'; ...
      'diag', @Rmem_diag, '(R)'; ...
      'denom', @Rmem_denom1, '(R, x)'; ...
      'denom_max0', @Rmem_denom_max0, '(R, x)'; ...
      'numer_pl_pcg_qs_ls', @Rmem_numer_pl_pcg_qs_ls, '(R, x1,x2)'; ...
      'denom_pl_pcg_qs_ls', @Rmem_denom_pl_pcg_qs_ls, '(R, x1,x2)'; ...
      'numer_denom_pl_pcg_qs_ls', @Rmem_numer_denom_pl_pcg_qs_ls, '(R, x1, x2)'; ...
      'cgrad_denom', @Rmem_cgrad_denom, '(R, x)'}; % not tested
  R = strum(R, meth);
else
  meth = {...
      'penal', @Rmem_penal, '(R, x)'; ...
      'cgrad', @Rmem_cgrad1_integ, '(R, x)'; ...
      'denom', @Rmem_denom1_integ, '(R, x)'; ...
      'cgrad_denom', @Rmem_cgrad_denom1_integ, '(R, x, numthreads)'; ...
      'denom_max0', @Rmem_denom_max0, '(R, x)'};
  R = strum(R, meth);
end

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
     if ~R.use_integ
	cgrad(:,ll) = Rmem_cgrad1(R, x);
     else
	cgrad(:,ll) = Rmem_cgrad1_integ(R, x);
     end
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

% $$$ old Rmem_denom_max0()
% $$$ %
% $$$ % Rmem_denom_max0()
% $$$ %
% $$$ function denom_max0 = Rmem_denom_max0(R)
% $$$ switch R.type_denom
% $$$ case {'matlab', 'aspire'}
% $$$ 	error 'not done due to R.wt'
% $$$ 	t = R.wt .* R.pot.wpot(R.pot, 0);
% $$$ 	if R.C_is_I
% $$$ 		denom_max0 = t;
% $$$ 	else
% $$$ 		t = reshape(t, [R.dim R.M]);
% $$$ 		if R.order == 1 % fix: order=2 denom?
% $$$ 			t = 2 * t; % "2" because (1,-1) differences
% $$$ 			denom_max0 = penalty_mex('diff1,back2', ...
% $$$ 					single(t), R.offsets);
% $$$ 		elseif R.order == 2 % (-1,2,1) so ?
% $$$ 			t = 4 * t; % 4 = |-1| + |2| + |-1|
% $$$ 			denom_max0 = penalty_mex('diff2,backA', ...
% $$$ 					single(t), R.offsets);
% $$$ 			warning 'todo: order=2 not tested'
% $$$ 		else
% $$$ 			error 'order not done for denom'
% $$$ 		end
% $$$ 		denom_max0 = R.denom_max0(R.mask);
% $$$ 	end
% $$$ otherwise
% $$$ 	error('Unknown type_denom: "%s"', R.type_denom)
% $$$ end

% new Rmem_denom_max0()
%
% Rmem_denom_max0()
%
function denom_max0 = Rmem_denom_max0(R,x)
if R.order==2
  error 'done for 1st order penalties only'
end
switch R.type_denom
case {'matlab', 'aspire'}
     denom_max0 =0;
     for mm=1:R.M
       % extra fwd proj for conservatively counting the neighbors
       d = penalty_mex_call(R.diff_str_forw, x, R.offsets(mm), 2);
       wpot0 = R.pot{mm}.wpot(R.pot{mm},0);
       d = 2 * R.beta(mm) * wpot0 * ( (R.kappa(R.mask).^2) .* ones(size(d)));
       d = penalty_mex('diff1,back2', single(d), R.offsets(mm));
       denom_max0 = denom_max0 +  d;
     end
otherwise
	error('Unknown type_denom: "%s"', R.type_denom)
end

%
% Rmem_numer_pl_pcg_qs_ls()
%
function numer = Rmem_numer_pl_pcg_qs_ls(R,x1,x2)
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
function denom = Rmem_denom_pl_pcg_qs_ls(R,x1,x2)
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
% function out = Rmem_numer_denom_pl_pcg_qs_ls(R,x1,x2)
function out = Rmem_numer_denom_pl_pcg_qs_ls(R,x1,x2)

if (prod(size(x1))==prod(size(col(R.mask==true)))) || ...
      (prod(size(x1))==prod(size(col(R.mask))))
  x_is_Cx = false;
else
  x_is_Cx = true;
  seglen = prod(size(x1))/R.M;
end

numer = 0;
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	if ~x_is_Cx
	  d1 = penalty_mex_call(R.diff_str_forw, x1, R.offsets(mm), 2); % 2 because a single column
	  d2 = penalty_mex_call(R.diff_str_forw, x2, R.offsets(mm), 2); % 2 because a single column
	else

	  d1 = reshape( x1( (mm-1)*seglen+1:mm*seglen ),size(R.kappa2) );
	  d2 = reshape( x2( (mm-1)*seglen+1:mm*seglen ),size(R.kappa2) );
	end
	d1_numer = R.beta(mm) * R.kappa2 .* pot.dpot(pot, d1);
	d1_denom = R.beta(mm) * R.kappa2 .* pot.wpot(pot, d1);
	numer = numer + sum(col(d1_numer .* d2));
	denom = denom + sum(col(d1_denom .* (d2.^2)));
end
out = [numer denom];

%
% Rmem_cgrad_denom
%
function out= Rmem_cgrad_denom(R, dummy, x)

if (prod(size(x))==prod(size(col(R.mask==true)))) || ...
      (prod(size(x))==prod(size(col(R.mask))))
  x_is_Cx = false;
else
  x_is_Cx = true;
  seglen = prod(size(x))/R.M;
end

cgrad = 0;
denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	if ~x_is_Cx
	  d = penalty_mex_call(R.diff_str_forw, single(x), R.offsets(mm), 2);
	else
	  d = reshape( x( (mm-1)*seglen+1:mm*seglen ),size(R.kappa2) );
	end
	d1 = R.beta(mm) * R.kappa2 .* pot.dpot(pot, d);
	cgrad = cgrad + penalty_mex_call(R.diff_str_back, single(d1), ...
		R.offsets(mm), 2); % 2 because a single column;
	
	d1 = 2 * R.beta(mm) * R.kappa2 .* pot.wpot(pot, d); % see denom code
	denom = denom + penalty_mex('diff1,back2', single(d1), R.offsets(mm));
end
out{1} = cgrad;
out{2} = denom;

%
% Rmem_denom1
%
function denom= Rmem_denom1(R, dummy, x)
siz = size(x);
flag_column = 0;
if size(x,1) == R.np % unmasked!!
	x = embed(x, R.mask);
	flag_column = 1;
end
x = reshape(x, prod(R.dim), []);

denom = 0;
for mm=1:R.M
	pot = R.pot{mm};
	d = penalty_mex_call(R.diff_str_forw, (single(x)), R.offsets(mm), ...
			     2);
	d1 = 2 * R.beta(mm) * R.kappa2 .* pot.wpot(pot, d); % see denom code
	denom = denom + penalty_mex('diff1,back2', single(d1), R.offsets(mm));
end

if flag_column
	denom = denom(R.mask(:), :);
end
denom = reshape(denom, siz);

%
% Rmem_dercurv()
% evaluate \dpoti and \wpoti
%
function [deriv, curv] = Rmem_dercurv(R, C1x)
error 'not done, ask jeff'
deriv = R.wt .* R.pot.wpot(R.pot, C1x) .* C1x; %ask:- cant we just do
                                               %dpot here?
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
% Rmem_offsets()
%
function offsets = Rmem_offsets(offsets, kappa);
if isempty(offsets)
	if ndims(kappa) == 2
		[nx ny] = size(kappa);
		offsets = [1 nx nx+1 nx-1];
	elseif ndims(kappa) == 3
		[nx ny nz] = size(kappa);
		offsets = [1 nx nx*ny];
	else
		error 'only 2D and 3D done'
	end
elseif streq(offsets, '3d:26') % all 26 neighbors (13 pairs)
	[nx ny nz] = size(kappa);
	offsets = [1 nx+[0 1 -1] nx*ny+col(outer_sum([-1:1],[-1:1]*nx))'];
end
offsets = int32(offsets);


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

%         _                    _
% _____/\/ %%%%%%%%%%%%%%%%%%%% \___     ___/|
%)_____  functions for use_integ ___)))))___ |
%      \/\_%%%%%%%%%%%%%%%%%%%%_/           \|
%

%
% Rmem_cgrad1_integ()
%
function cgrad = Rmem_cgrad1_integ(R, dummy, x)

siz = size(x);
flag_column = 0;
if size(x,1) == R.np % unmasked!!
	x = embed(x, R.mask);
	flag_column = 1;
end
x = reshape(x, prod(R.dim), []);

temp = R.pot_arg{1};
pot_params = zeros([size(temp,2)-1 1]);
for ii=1:size(temp,2)-1
  pot_params(ii) = temp{ii+1}(1);
end
cgrad = penalty_mex(sprintf('cgrad%d',R.order), ...
  		    single(x), single(R.beta), single(R.kappa),  ...
		    int32(R.offsets), temp{1}, single(pot_params) );

if flag_column
	cgrad = cgrad(R.mask(:), :);
end
cgrad = reshape(cgrad, siz);

%
% Rmem_denom1
%
function denom= Rmem_denom1_integ(R, dummy, x)
siz = size(x);
flag_column = 0;
if size(x,1) == R.np % unmasked!!
	x = embed(x, R.mask);
	flag_column = 1;
end
x = reshape(x, prod(R.dim), []);

temp = R.pot_arg{1};
pot_params = zeros([size(temp,2)-1 1]);
for ii=1:size(temp,2)-1
  pot_params(ii) = temp{ii+1}(1);
end
denom = penalty_mex(sprintf('denom%d',R.order), ...
  		    single(x), single(R.beta), single(R.kappa),  ...
		    int32(R.offsets), temp{1}, single(pot_params) );

if flag_column
	denom = denom(R.mask(:), :);
end
denom = reshape(denom, siz);

%
% Rmem_cgrad_denom1_integ
%
function cgrad_denom=Rmem_cgrad_denom1_integ(R, dummy, x, numthreads)
siz = size(x);
sizo = siz;
sizo(1) = 2*sizo(1);

flag_column = 0;
if size(x,1) == R.np % unmasked!!
	x = embed(x, R.mask);
	flag_column = 1;
end
x = reshape(x, prod(R.dim), []);

temp = R.pot_arg{1};
pot_params = zeros([size(temp,2)-1 1]);
for ii=1:size(temp,2)-1
  pot_params(ii) = temp{ii+1}(1);
end
cgrad_denom = penalty_mex(sprintf('cgrad_denom%d',R.order), ...
  		    single(x), single(R.beta), single(R.kappa),  ...
		    int32(R.offsets), temp{1}, single(pot_params), numthreads );

if flag_column
        temp1 = cgrad_denom(1:size(cgrad_denom,1)/2);
        temp2 = cgrad_denom(size(cgrad_denom,1)/2+1:end);
	cgrad_denom = [temp1(R.mask(:), :) temp2(R.mask(:), :)];
end
cgrad_denom = reshape(cgrad_denom, sizo);


