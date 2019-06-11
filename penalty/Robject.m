 function R = Robject(kappa, varargin)
%function R = Robject(kappa, [options])
%|
%| Build roughness penalty regularization "object" based on C = Cdiff() object,
%| for regularized solutions to inverse problems.
%|
%| General form of nonquadratic penalty function:
%|	R(x) = \sumk w_k \pot([Cx]_k), where [Cx]_k = \sum_j c_{kj} x_j.
%|	
%| For quadratic case, \potk(t) = t^2/2 in which case
%|	R(x) = x' C' W C x / 2, where W depends on beta and edge_type.
%|
%| Penalty gradient is C' D C x, where D = diag{\wpot_k([Cx]_k)}.
%|
%| in
%|	kappa	[nx,ny[,nz]]	kappa array, or logical support mask
%|
%| options
%|	edge_type '?'		how to handle mask edge conditions (see Cdiff)
%|		'none'		no roughness penalty (NOT DONE)
%|		'tight'		only penalize within-mask differences (default)
%|		'leak'		penalize between mask and neighbors
%|				(its primary use for consistency with ASPIRE)
%|	'order', ?		1st-order or 2nd-order differences (see Cdiff)
%|	'offsets', [?]		offsets to neighboring pixels
%|					(see Cdiff for the defaults)
%|				use '3d:26' to penalize all 13 pairs of nbrs
%|				use '0' for C = I (identity matrix)
%|	'beta', ?		global regularization parameter
%|					default: 2^0
%|	'delta', ?		potential parameter, see potential_func()
%|					or {delta, param}.  default: inf
%|	'potential', '?'	e.g., 'huber', see potential_func()
%|				default: 'quad' for quadratic regularization.
%|	'type_denom', '?'	type of "denominator"
%|					(for quadratic surrogates like SPS)
%|		'matlab'	denominator for SPS
%|		'aspire'	denominator for SPS that matches aspire
%|		'none'		no denominator (default)
%|				(because R.E and R.denom needed only for SPS)
%|	'distance_power', ?	1 classical (default), 2 possibly improved
%|					see penalty_mex()
%|	'user_wt', [?]		User-provided array of penalty weight values
%|					size: [nx,ny[,nz],#offsets]
%|				of dimension [size(mask) length(offsets)].
%|				These are .* the usual wt values for edge_type.
%|	'mask'			Override default: mask = (kappa ~= 0)
%|
%| out
%|	R structure has the following "methods" 
%|	R.penal(R, x)	evaluates R(x)
%|	R.cgrad(R, x)	evaluates \cgrad R(x) (column gradient)
%|	R.denom(R, x)	evaluates denominator for separable surrogate
%|	[pderiv pcurv] = feval(R.dercurv, R, C1*x) derivatives and curvatures
%|			for non-separable parabola surrogates
%|	R.diag(R)	diagonal of Hessian of R (at x=0), for preconditioners.
%|	R.C1		differencing matrix, with entries 1 and -1,
%|			almost always should be used in conjunction with R.wt
%|
%| Typical use:	mask = true([128 128]); % or something more conformal
%|		R = Robject(mask, 'beta', 2^7);
%|
%| Copyright 2004-11-14, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(kappa, 'test'), run_mfile_local 'Robject_test', return, end

% option defaults
R.potential = 'quad';
R.beta = 2^0;
R.delta = inf;
R.edge_type = 'tight';
R.type_denom = 'none';
R.distance_power = 1;
R.order = 1; % 1st-order differences, used in call to Cdiff
R.user_wt = []; % user-provided wt values
R.mask = [];
R.offsets = [];

% parse name/value option pairs
R = vararg_pair(R, varargin);

% dimensions, and default offsets
if isempty(R.offsets)
	if ndims(kappa) == 2
		[nx ny] = size(kappa);
		R.offsets = [1 nx nx+1 nx-1];
	elseif ndims(kappa) == 3
		[nx ny nz] = size(kappa);
		R.offsets = [1 nx nx*ny];
	else
		error 'only 2D and 3D done'
	end
elseif streq(R.offsets, '3d:26') % all 26 neighbors (13 pairs)
	[nx ny nz] = size(kappa);
	R.offsets = [1 nx+[0 1 -1] nx*ny+col(outer_sum([-1:1],[-1:1]*nx))'];
end

if ~iscell(R.delta)
	R.delta = {R.delta};
end

R.offsets = int32(R.offsets);

if streq(R.edge_type, 'none')
	error 'todo: unpenalized'
end

if R.beta < 0, warn('Negative beta? This is probably wrong!'), end

R.isquad = streq(R.potential, 'quad');

if isempty(R.mask)
	R.mask = kappa ~= 0; % default is to infer from kappas
end
if ~islogical(R.mask), error 'mask must be logical', end

%
% build plain sparse differencing object "C1", containing 1 and -1 values
% (or identity matrix)
%
R.C_is_I = any(R.offsets == 0);
if R.C_is_I
	if any(R.offsets), error 'identity is offsets=0', end
	R.C1 = diag_sp(ones(sum(R.mask(:)),1)); % identity matrix
	R.wt = R.beta;

else
	R.C1 = Cdiff(R.mask, 'offsets', R.offsets, ...
			'edge_type', 'none', 'order', R.order);

	% wk factors
	R.wt_string = sprintf('wk,%s,%d', R.edge_type, R.order);
	R.wt = penalty_mex(R.wt_string, single(kappa), R.offsets, R.distance_power);
	R.wt = R.beta * single(R.wt(:)); % absorb beta into wk factors
end

% incorporate user provided wk values
if ~isempty(R.user_wt)
	R.wt = single(R.wt .* R.user_wt(:));
	R.user_wt = []; % save memory
end

% trick: for quadratic case, provide a R.C object where C = sqrt(W) * C1
if R.isquad
	R.C = R.C1;
	R.C.cascade_after = diag_sp(sqrt(R.wt));
end

%
% desired potential function
%
R.pot = potential_func(R.potential, R.delta{:});

%
% functions
%
R.dercurv = @Robject_dercurv;
R.handle_denom = @Robject_denom;
R.handle_diag = @Robject_diag;
% trick: the following form of R.penal allows multiple realizations of x
R.penal = @(R, x) ... % 2014-04-27 added 'double' to better match new version
	sum(repmat(R.wt, ncol(x)) .* R.pot.potk(R.pot, R.C1 * x), 'double');
R.cgrad = @(R, x) R.C1' * (diag_sp(R.wt) * R.pot.dpot(R.pot, R.C1 * x));
R.diag = @(R, x) R.handle_diag(R);
R.denom = @(R, x) R.handle_denom(R, x);

%
% precompute denominator for separable quadratic case if needed
%
if ~streq(R.type_denom, 'none') && R.isquad
	if streq(R.type_denom, 'matlab') || streq(R.type_denom, 'aspire')
		t = R.wt .* R.pot.wpot(R.pot, 0);
		if R.C_is_I
			R.denom_max0 = t;
		else
			t = reshape(t, [size(R.mask) length(R.offsets)]);
			if R.order == 1 % fix: order=2 denom?
				t = 2 * t; % "2" because (1,-1) differences
				R.denom_max0 = penalty_mex('diff1,back2', ...
						single(t), R.offsets);
			elseif R.order == 2 % (-1,2,1) so ?
				t = 4 * t; % 4 = |-1| + |2| + |-1|
				R.denom_max0 = penalty_mex('diff2,backA', ...
						single(t), R.offsets);
			else
				error 'order not done for denom'
			end
			R.denom_max0 = single(R.denom_max0(R.mask));
		end
	else
		error(['Unknown type_denom: ' R.type_denom])
	end
end


%
% Robject_denom()
% jth penalty separable surrogate curvature is d_j = \sumk |\ckj| \ck \wpotk
% where \ck = \sumj |\ckj|.
% Here, ck = 2 since there is a +1 and a -1 per row of C1 (for 1st-order) 
% Also, |\ckj| = |\ckj|^2 since \ckj = +/- 1, so we can use 'diff1,back2'
%
function denom = Robject_denom(R, x)
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
		t = reshape(t, [size(R.mask) length(R.offsets)]);
		if R.order == 1
			t = 2 * t; % because (1,-1)
			denom = penalty_mex('diff1,back2', t, R.offsets);
		elseif R.order == 2 % (-1,2,1) so ck = 4
			t = 4 * t; % 4 = |-1| + |2| + |-1|
			denom = penalty_mex('diff2,backA', ...
					single(t), R.offsets);
		else
			error 'order not done'
		end
		denom = single(denom(R.mask));
	end
end


%
% evaluate diagonal of Hessian of (quadratic surrogate for) R at x=0
% \sum_k w_k |c_kj|^2 \wpot(0)
%
function rjj = Robject_diag(R)
t = single(R.wt .* R.pot.wpot(R.pot, 0));
t = reshape(t, [size(R.mask) length(R.offsets)]);
if R.C_is_I
	error 'not done, ask jeff'
else
	rjj = penalty_mex('diff1,back2', t, R.offsets);
	rjj = double(rjj(R.mask));
end


%
% evaluate \dpoti and \wpoti
%
function [deriv, curv] = Robject_dercurv(R, C1x)
deriv = R.wt .* R.pot.wpot(R.pot, C1x) .* C1x; 
curv = R.wt .* R.pot.wpot(R.pot, C1x);


%
% Compute both cgrad and denom (of separable surrogate) efficiently.
% Unused for now, but could be used if anonymous functions are too inefficient
% since both R.cgrad and R.denom use C*x so there is redundancy.
% What we really need is an anonymous function that has two output arguments.
% No, the feval with a function_handle will suffice!
%
%function [cgrad, denom] = Robject_cgrad_denom(R, x)
%if R.isquad
%	cgrad = R.C' * (R.C * x);
%	denom = R.denom;
%else
%	Cx = R.C * x;
%	wx = R.wt .* R.pot.wpot(R.pot, Cx);
%	cgrad = R.C' * (wx .* Cx);
%	denom = R.E * wx;
%end

