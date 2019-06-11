  function wt = Rweights(kappa, offsets, varargin)
%|function wt = Rweights(kappa, offsets, [options])
%|
%| Build "weights" for roughness penalty for regularized methods.
%| Intended to be used internally by roughness penalty object Reg1.
%|
%| General form of roughness penalty:
%|	R(x) = sum_{m=1}^M sum_n w[n;m] potential( [C_m x]_n )
%| where M is the number of neighbors (offsets), and n=1,...,N=numel(x)
%| and C_m is a (square) differencing matrix in the mth direction off_m.
%|
%| General form of weights:
%|	w[n;m] = beta_m / || off_m ||^distance_power kappa2[n;m] user_wt[n;m]
%| where form of kappa2[n;m] depends on 'edge_type' as follows:
%|	kappa^2[n] for 'simple' case,
%|	kappa[n] kappa[n-off_m] for 'tight', order = 1
%|	kappa[n] sqrt(kappa[n-off_m] * kappa[n+off_m]) for 'tight', order = 2
%|	kappa[n] kappa[n-off_m] for 'leak' case, order = 1,
%|			unless either is zero, in which case square the other.
%|			for order=2, square the maximum of all three if needed.
%|	kappa[n] kappa_extend[n-off_m] for 'aspire2' case, order = 1
%|
%| Although the 'tight' case seems preferable to avoid any influence
%| of pixels outside of the support, the 'simple' case has the advantage
%| that it can be precomputed easily with only N storage, instead of M*N.
%|
%| in
%|	kappa	[(N)]		kappa array, or logical support mask
%|	offsets	[M]		offset(s) to neighboring pixels
%|
%| options
%|	'type_wt'		what type of object to return:
%|		'array'		[(N) M] array of w[n;m] values
%|		'fly'		strum object for computing w[:,m] on-the-fly
%|		'pre'		strum object for precomputed w[:,m] (default)
%|	'edge_type'		how to handle mask edge conditions
%|		'simple'	kappa^2 (almost "tight" but saves memory)
%|		'tight'		only penalize within-mask differences (default)
%|		'leak'		penalize between mask and neighbors
%|					(mostly for consistency with ASPIRE)
%|		'aspire2'	aspire 2d version of 'leak' (order=1 only)
%|					(for jf testing only - not recommended)
%|	'beta', [1] or [M]	regularization parameter(s) (default: 2^0)
%|	'order' {1|2}		differencing order (only for some cases)
%|	'distance_power', {0:1:2}  1 classical (default), 2 possibly improved
%|					use 0 if user_wt has its effect already
%|	'user_wt', [(N) M]	User-provided array of penalty weight values
%|				of dimension [size(kappa) length(offsets)].
%|	'use_mex' {0|1}		use penalty_mex? (default: 0)
%|
%| out
%|	wt [(N) M]		w[n;m] values needed in regularizer
%|	or a strum object with the following public methods:
%|		wt.col(m)	returns w[:,m]
%|
%| Copyright 2006-12-2, Jeff Fessler, University of Michigan

if nargin == 1 && streq(kappa, 'test'), Rweights_self_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% arguments
arg.offsets = offsets;

% option defaults
arg.type_wt = 'pre'; % precompute [(N) M] array, but access via strum
arg.use_mex = false;
arg.beta = 2^0;
arg.edge_type = 'tight';
arg.distance_power = 1;
arg.order = 1; % 1st-order differences
arg.user_wt = []; % user-provided wt values

% parse name/value option pairs
arg = vararg_pair(arg, varargin);

% dimensions
arg.isize = size(kappa);

% beta values adjusted for distances
arg.beta_dis = arg.beta(:) ...
	./ penalty_distance(arg.offsets(:), arg.isize) .^ arg.distance_power;

arg.displace = penalty_displace(arg.offsets, arg.isize); % [M *N]

kappa = single(kappa); % no logical

% 'aspire2' option is only for testing the special case of 2d.
if streq(arg.edge_type, 'aspire2')
	nx = arg.isize(1);
	if arg.order ~= 1 || ndims(kappa) ~= 2 || arg.distance_power ~= 2 ...
		|| arg.use_mex || ~isequal(arg.offsets, [1 nx nx+1 nx-1])
		error 'aspire2 needs 2d, order=1, offset=[1 nx nx+1 nx-1]'
	end
	if any(col(kappa([1 end],:))) || any(col(kappa(:,[1 end])))
		error 'aspire2 requires support exclude one pixel border'
	end
end

switch arg.type_wt
case 'array'
	wt = Rweights_array(arg, kappa);

case 'fly'
	wt = Rweights_strum_fly(arg, kappa);

case 'pre'
	wt = Rweights_strum_pre(arg, kappa);

otherwise
	fail('unknown type %s', arg.type_wt)
end

end % Rweights


% Rweights_strum_fly()
function wt = Rweights_strum_fly(arg, kappa)
if streq(arg.edge_type, 'simple')
	arg.kappa2 = kappa(:).^2;
	fun = @Rweights_strum_fly_col_simple;
else
	arg.kappa = kappa;
	if arg.use_mex
		fun = @Rweights_strum_fly_col_mex;
	else
		fun = @Rweights_strum_fly_col_array;
	end
end
if ~isempty(arg.user_wt)
	arg.fun_no_user_wt = fun;
	fun = @Rweights_strum_mul_user_wt;
end
wt = {...
	'col', fun, '(l): returns w(:,l)'; ...
	'show', @Rweights_strum_show, '([im opts])'
	};
wt = strum(arg, wt);
end % Rweights_strum_fly()


% Rweights_strum_show()
function out = Rweights_strum_show(arg, varargin)
Nd = arg.isize;
Ns = prod(Nd);
MM = length(arg.offsets);
tmp = zeros(Ns, MM);
for mm=1:MM
	tmp(:,mm) = arg.col(mm);
end
tmp = reshapee(tmp, [Nd MM]);
im(tmp, varargin{:})
out = [];
end % Rweights_strum_show()


% Rweights_strum_pre()
function wt = Rweights_strum_pre(arg, kappa)
arg.w = Rweights_array(arg, kappa); % user_wt gets built in here!
fun = @ (arg,mm) arg.w(:,mm);
tmp = {...
	'col', fun, '(l): returns w(:,l)';
	'show', @Rweights_strum_show, '([im opts])'
	};
wt = strum(arg, tmp);
end % Rweights_strum_pre()


% Rweights_strum_mul_user_wt()
function wt = Rweights_strum_mul_user_wt(arg, mm)
wt = arg.fun_no_user_wt(arg, mm);
tmp = stackpick(arg.user_wt, mm);
wt = wt .* col(tmp);
end % Rweights_strum_mul_user_wt()


% Rweights_strum_fly_col_array()
function wt = Rweights_strum_fly_col_array(arg, mm)
%if ~isvar('mm'), mm = 1:length(arg.offsets); end % default: return all
wt = arg.beta_dis(mm) * Rweights_kappa2_mat(arg.kappa, arg.offsets(mm), ...
	arg.displace(mm,:), arg.edge_type, arg.order); % [N* 1]
end % Rweights_strum_fly_col_array()


% Rweights_strum_fly_col_mex()
function wt = Rweights_strum_fly_col_mex(arg, mm)
wt = arg.beta_dis(mm) * Rweights_kappa2_mex(arg.kappa, arg.offsets(mm), ...
		arg.edge_type, arg.order, 0);
wt = wt(:);
end % Rweights_strum_fly_col_mex()


% Rweights_strum_fly_col_simple()
function wt = Rweights_strum_fly_col_simple(arg, mm)
wt = arg.beta_dis(mm) * arg.kappa2;
end % Rweights_strum_fly_col_simple()


% Rweights_array()
function wt = Rweights_array(arg, kappa)
MM = length(arg.offsets);

if arg.use_mex
	try
		wt = Rweights_kappa2_mex(kappa, arg.offsets, arg.edge_type, ...
			arg.order, 0);
		wt = reshapee(wt, [], MM);
		for mm=1:MM
			wt(:,mm) = arg.beta_dis(mm) * wt(:,mm);
		end

		if ~isempty(arg.user_wt) % incorporate user provided weights
			tmp = reshape(arg.user_wt, size(wt));
			wt = single(wt .* tmp);
		end
	return
	catch
		warn 'mex unavailable despite request; resorting to matlab'
	end
end

wt = zeros(prod(arg.isize), MM);
for mm=1:MM
	wt(:,mm) = arg.beta_dis(mm) * Rweights_kappa2_mat(kappa, ...
		arg.offsets(mm), arg.displace(mm,:), arg.edge_type, arg.order);
end

if ~isempty(arg.user_wt) % incorporate user provided weights
	tmp = reshape(arg.user_wt, size(wt));
	wt = single(wt .* tmp);
end

end % Rweights_array()


% Rweights_kappa2_mex()
function wt = Rweights_kappa2_mex(kappa, offsets, edge_type, order, distance_power)
wt_string = sprintf('wk,%s,%d', edge_type, order);
wt = penalty_mex(wt_string, single(kappa), int32(offsets), distance_power);
end % Rweights_kappa2_mex()


% Rweights_kappa2_mat()
% output is [N* L]
function kappa2 = Rweights_kappa2_mat(kappa, offset, displace, edge_type, order)

pn = jf_protected_names;

Ns = numel(kappa);
Nd = size(kappa);

switch edge_type
case 'simple'
	kappa2 = kappa(:).^2;

case 'leak'
	kappa2 = zeros(size(kappa(:)));
	if order == 1
		ii = (1+max(offset,0)):(Ns+min(offset,0));
		kappa0 = kappa(ii) + kappa(ii-offset) .* ~kappa(ii);
		kappa1 = kappa(ii-offset) + kappa(ii) .* ~kappa(ii-offset);
		kappa2(ii) = kappa0 .* kappa1;

		if 1 % handle edge effects properly!
			ic = pn.ind2sub(Nd, ii);
			in = pn.ind2sub(Nd, ii-offset);
			good = (ic - in) == repmat(displace, [length(ii) 1]);
			kappa2(ii) = kappa2(ii) .* all(good, 2);
		end

	elseif order == 2
		ii = [(1+abs(offset)):(Ns-abs(offset))]';
		kaps = [kappa(ii) kappa(ii-offset) kappa(ii+offset)];
		tmp = kaps(:,1) .* sqrt(kaps(:,2) .* kaps(:,3));
		bad = tmp == 0;
		tmp(bad) = max(kaps(bad,:), [], 2).^2;
		kappa2(ii) = tmp;

		if 1 % handle edge effects properly!
			ic = pn.ind2sub(Nd, ii);
			in = pn.ind2sub(Nd, ii-offset);
			ip = pn.ind2sub(Nd, ii+offset);
			tmp = repmat(displace, [length(ii) 1]);
			good = ((ic - in) == tmp) & ((ip - ic) == tmp);
			kappa2(ii) = kappa2(ii) .* all(good, 2);
		end
	else
		fail('bad order %d', order)
	end

case 'tight'
	kappa2 = zeros(size(kappa(:)));
	if order == 1
		ii = (1+max(offset,0)):(Ns+min(offset,0));
		kappa2(ii) = kappa(ii) .* kappa(ii-offset);

		if 1 % handle edge effects properly!
			ic = pn.ind2sub(Nd, ii);
			in = pn.ind2sub(Nd, ii-offset);
			good = (ic - in) == repmat(displace, [length(ii) 1]);
			kappa2(ii) = kappa2(ii) .* all(good, 2);
		end

	elseif order == 2
		ii = (1+abs(offset)):(Ns-abs(offset));
		kappa2(ii) = kappa(ii) ...
			.* sqrt(kappa(ii-offset) .* kappa(ii+offset));

		if 1 % handle edge effects properly!
			ic = pn.ind2sub(Nd, ii);
			in = pn.ind2sub(Nd, ii-offset);
			ip = pn.ind2sub(Nd, ii+offset);
			tmp = repmat(displace, [length(ii) 1]);
			good = ((ic - in) == tmp) & ((ip - ic) == tmp);
			kappa2(ii) = kappa2(ii) .* all(good, 2);
		end

	else
		fail('bad order %d', order)
	end

case 'aspire2' % match rp_beta_fix_edges() in rp,beta.c
	if order ~= 1 || ndims(kappa) ~= 2, error 'bug', end
	kappa = Rweights_kappa_expand(kappa);
	kappa2 = Rweights_kappa2_mat(kappa, offset, displace, 'tight', order); % trick

otherwise
	fail('unknown edge_type "%s"', edge_type)
end
end % Rweights_kappa2_mat()


% Rweights_kappa_expand()
% match rp_beta_fix_edges() in rp,beta.c
function kappa = Rweights_kappa_expand(kappa)

if any(col(kappa([1 end],:))) || any(col(kappa(:,[1 end])))
	error 'kappa_expand requires support one pixel from edge'
end
mask = kappa ~= 0;
[nx ny] = size(mask);
	for iy=2:(ny-1)
	 for ix=2:(nx-1)
		if ~mask(ix,iy), continue, end
		for dy = -1:1
		 for dx = -1:1
			if ~mask(ix+dx, iy+dy)
				kappa(ix+dx, iy+dy) = kappa(ix,iy);
			end
 		 end
		end
	 end
	end
end % Rweights_kappa_expand()


% mask_border_check()
% check if mask has nonzero values along its edges.
% this seems unnecessary because of how the code has been written.
function mask_border_check(mask)
for id=1:ndims(mask)
	tmp = shiftdim(mask, id-1);
	tmp = reshape(tmp, size(tmp,1), []);
	if any(tmp(1,:)) || any(tmp(end,:))
		warn('mask is nonzero at boundaries along dim=%d', id)
	end
end
end % mask_border_check()



% Rweights_self_test
function Rweights_self_test
if 0
	ig = image_geom('nx', 10, 'ny', 8, 'nz', 6, 'dx', 1, 'dz', 2);
	offsets = [1 ig.nx -ig.nx ig.nx+1 ig.nx-1 ig.nx*ig.ny ig.nx*ig.ny+1];
else
	ig = image_geom('nx', 10, 'ny', 8, 'dx', 1);
	offsets = [1 ig.nx+[0 1 -1]]; % -ig.nx+1];
end
ig.mask = ig.circ > 0;
kappa = ig.embed(1+[1:ig.np]'/ig.np);

rng(0)
user_wt = 5+randn([ig.dim length(offsets)]);
uwa = {'user_wt', user_wt};

% test both matlab and mex versions
edge_types = {'tight', 'leak', 'simple'};
im plc 2 4
for use_wt = [1 0] % first with user_wt, then without for figures
	if ~use_wt
		uwa = {}; % conventional
	end
for order=1:2
for ie=1:length(edge_types)
	edge_type = edge_types{ie};
	printm('use_wt=%d order=%d edge_type=%s', use_wt, order, edge_type)
	arg = {kappa, offsets, 'order', order, ...
		'type_wt', 'array', 'edge_type', edge_type, uwa{:}};
	wt = Rweights(arg{:});

	if ~streq(edge_type, 'simple') && has_mex_jf
		wx = Rweights(arg{:}, 'use_mex', 1);
		equivs(wt, wx)
	end
	im(ie+4*(order-1), ig.shape(wt))
	xlabelf('order=%d type=%s', order, edge_type)

	wp = Rweights(arg{:}, 'type_wt', 'pre');
	ws = Rweights(arg{:}, 'type_wt', 'fly');
	wy = Rweights(arg{:}, 'type_wt', 'fly', 'use_mex', 1);

	for mm=1:ncol(wt)
		jf_equal(wp.col(mm), wt(:,mm))
		equivs(ws.col(mm), wt(:,mm))

		if has_mex_jf && max_percent_diff(wt(:,mm), wy.col(mm)) > 2e-5
			im plc 1 3
			im(1, ig.shape(wt(:,mm)), 'wt'), cbar
			im(2, ig.shape(wy.col(mm)), 'wy'), cbar
			im(3, ig.shape(wy.col(mm) - wt(:,mm)), 'wy'), cbar
			error 'mex vs mat bug'
		end
	end
end % ie
end % order
end % use_wt
im subplot 8
ws.show('test show')

% test aspire2 special case
kappa(:, [1 end]) = 0; kappa([1 end], :) = 0; % no border support for aspire
arg = {kappa, [1 ig.nx+[0 1 -1]], 'order', 1, ...
	'type_wt', 'array', 'edge_type', 'aspire2', 'distance_power', 2};
wt = Rweights(arg{:});
im(4, ig.shape(wt)), xlabel 'aspire2'
ws = Rweights(arg{:}, 'type_wt', 'fly');
wp = Rweights(arg{:}, 'type_wt', 'pre');
for mm=1:ncol(wt)
	jf_equal(ws.col(mm), wt(:,mm))
	jf_equal(wp.col(mm), wt(:,mm))
end

end % Rweights_self_test
