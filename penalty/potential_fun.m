 function pot = potential_fun(ptype, delta, param, varargin)
%function pot = potential_fun(ptype, delta, param, varargin)
%|
%| Define roughness penalty potential functions (as strum).
%|
%| The penalty will have the form
%|	R(x) = sum_k w_k * potential_k([Cx]_k, delta_k)
%| where w_k is provided elsewhere, not here!
%|
%| in
%|	ptype		quad broken huber hyper2 hyper3 cauchy qgg2 genhub
%|			lange1 lange3 (fair) geman&mcclure gf1 li98cfs
%|			Recommended: 'hyper3'
%|			To list all possible choices, use:
%|			potential_fun('list') for "smooth" options
%|			potential_fun('list1') for "non-smooth" options
%|	delta		scalar, or image-sized array;
%|			"cutoff" parameter for edge-preserving regularization
%|	param		optional additional parameter(s) for some choices:
%|				'gf1' (generalized Fair) : [a b]
%|				'qgg2' : q
%|				'genhub' & 'stevenson94dpr' : [p q]
%|				'table1' : [dz, dpot([1:K] * dz)]
%|
%| option
%|	'dummy'		include extra dummy argument for backward compatibility
%|			e.g., pot.potk([], C*x). default: 0
%|
%| out
%|	pot		strum object, with data: delta and param
%|	methods:
%|		pot.potk(C*x)	potential function value
%|		pot.wpot(C*x)	potential 'weights' (aka half-quad. curvatures)
%|		pot.dpot(C*x)	potential derivative
%| 		pot.shrink(b, reg)	proximal (shrinkage) operation:
%|					argmin_z 1/2 |z - b|^2 + reg * pot(z)
%|		pot.plot()	plot all of the above functions
%|
%| trick: for ptype 'gf1-fit', the param argument should be:
%|	{'name of potential to fit', points, param}
%| and this function returns the [a b] parameters needed for a subsequent
%| call with ptype 'gf1'
%|
%| Copyright 2004-5-18, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(ptype, 'test') ir_potential_fun_test, return, end
if streq(ptype, 'test_lpp') ir_potential_fun_test_lpp, return, end

if streq(ptype, 'list') % return list of all "smooth" choices
	pot = {'quad', 'huber', 'huber2', ...
	'hyper2', 'hyper3', 'cauchy', ...
	'qgg2', 'genhub', 'lange1', 'lange3', 'geman&mcclure', 'gf1', ...
	'table0', 'table1', 'li98cfs'};
%	'stevenson94dpr' omit this one; use 'genhub' instead
return
end

if streq(ptype, 'list1') % return list of 'l1' type choices (mostly non-smooth)
	pot = {'l0', 'l1', 'lpp', 'tav', 'broken', 'fair-l1'};
return
end

if ~isvar('delta'), delta = []; end
if ~isvar('param'), param = []; end

arg.dummy = false;
arg.scale = 1;
arg = vararg_pair(arg, varargin);

if streq(ptype, 'gf1-fit') % trick: just return parameters for this case
	pot = ir_potential_fun_gf1_fit(param{:});
return
end

[pot.type, pot.delta, pot.param, potk, wpot, dpot, shrink] = ...
	ir_potential_fun_parse(ptype, delta(:), param(:));

if arg.scale ~= 1
	potk = @(pot, z) arg.scale * potk(pot, z);
	dpot = @(pot, z) arg.scale * dpot(pot, z);
	wpot = @(pot, z) arg.scale * wpot(pot, z);
end

if arg.dummy
	potk = @(pot, dum, z) potk(pot, z);
	dpot = @(pot, dum, z) dpot(pot, z);
	wpot = @(pot, dum, z) wpot(pot, z);
end

meth = {'potk', potk, 'wpot', wpot, 'dpot', dpot, 'shrink', shrink, ...
	'plot', @ir_potential_fun_plot};
pot = strum(pot, meth);

end % potential_fun()


% ir_potential_fun_gf1_fit()
% trick: gf1-fit means match gf1 to named potential at given points [0, s1, s2]
% where the given point values are s = z / delta, i.e., in "delta" units
% param: {'name of potential to fit', points, param}
function ab = ir_potential_fun_gf1_fit(ptype, sv, param)
sv = sv(:);
pot = potential_fun(ptype, 1, param); % trick: delta=1 here
pt = pot.wpot(sv);

% ab = [1 1; -pt']' \ ((pt-1) ./ sv);

s1 = sv(1);
s2 = sv(2);
w1 = pt(1);
w2 = pt(2);
ab(1) = (w2 * s2 * (1 - w1) - w1 * s1 * (1 - w2)) ...
		/ ((w1 - w2) * s1 * s2);
ab(2) = (s2 * (1 - w1) - s1 * (1 - w2)) ...
		/ ((w1 - w2) * s1 * s2);
end % ir_potential_fun_gf1_fit()


% ir_potential_fun_parse()
function [ptype, delta, param, potk, wpot, dpot, shrink] = ...
	ir_potential_fun_parse(ptype, delta, param)

dpot = [];
shrink = [];

% trick: huber2 is just a huber with delta / 2
% so that weighting function drops to 1/2 at delta, like hyper3 etc.
if streq(ptype, 'huber2')
	ptype = 'huber';
	delta = delta / 2;
end

% trick: hyper3 is just a hyperbola with delta scaled by sqrt(3)
% to approximately "match" the properties of 'cauchy' (old erroneous 'hyper')
if streq(ptype, 'hyper3')
	ptype = 'hyper2';
	delta = delta / sqrt(3);
end

if streq(ptype, 'gf1') && param(1) == 0
	if param(2) ~= 1
		fail 'only b=1 makes sense for gf1 with a=0'
	end
	ptype = 'lange3'; param = []; % trick: gf1 with a=0 and b=1 is lange3
end

switch ptype

% quadratic potential function
case 'quad'
	potk = @(pot, z) (abs(z).^2) / 2;
	wpot = @(pot, z) ones(size(z), class(z));
	dpot = @(pot, z) z;
	shrink = @(pot, b, reg) b ./ (1 + reg);

% broken parabola
case 'broken'
	potk = @(pot, z) min(z.^2, pot.delta.^2)/2;
	wpot = @(pot, z) ir_tonumeric(abs(z) < pot.delta, z);
	dpot = @(pot, z) z .* (abs(z) < pot.delta);
	shrink = @(pot, b, reg) ir_broken_shrink(b, reg, pot.delta);

% l0 (hard threshold)
case 'l0'
	potk = @(pot, z) ir_tonumeric(z ~= 0, z);
	wpot = @(pot, z) nan(size(z), class(z));
	dpot = @(pot, z) nan(size(z), class(z));
	shrink = @(pot, b, reg) b .* (abs(b) > sqrt(2*reg));

% absolute value (i.e., l1, for soft threshold)
case {'abs', 'l1'}
	potk = @(pot, z) abs(z);
	wpot = @(pot, z) 1 ./ abs(z); % invalid at z=0
	dpot = @(pot, z) sign(z);
	shrink = @(pot, b, reg) sign(b) .* max(abs(b) - reg, 0);

% l_p^p aka generalized gaussian (gg)
case 'lpp'
	potk = @(pot, z) abs(z).^param;
	dpot = @(pot, z) param * sign(z) .* abs(z).^(param-1); % invalid at z=0
	wpot = @(pot, z) param .* abs(z).^(param-2); % invalid at z=0
	switch param
	case 0
		potk = @(pot, z) ir_tonumeric(z ~= 0, z);
		dpot = @(pot, z) zeros(size(z), class(z)); % meaningless
		shrink = @(pot, b, reg) b .* (abs(b) > sqrt(2*reg));
	case 1
		potk = @(pot, z) abs(z);
		dpot = @(pot, z) sign(z); % invalid at t=0
		shrink = @(pot, b, reg) sign(b) .* max(abs(b) - reg, 0);
	case 0.5 % p=1/2
		shrink = @(pot, b, reg) ir_lpp_1_2_shrink(b, reg);
	case 4/3. % p=4/3
		shrink = @(pot, b, reg) ir_lpp_4_3_shrink(b, reg);
	case 1.5 % p=3/2
		shrink = @(pot, b, reg) sign(b) .* ...
			( sqrt((3/4*reg).^2 + abs(b)) - 3/4*reg ).^2;
		% from eqn. (4.5) in chaux:07:avf
	%	z + 9 * reg.^2 * sign(z) .* (1 - sqrt(1 + 16*abs(z)/(9*reg^2))) / 8;

	otherwise
		warn('no shrink for p=%g', param)
	end

% truncated absolute value
case 'tav'
	potk = @(pot, z) min(abs(z), pot.delta);
	wpot = @(pot, z) (1 ./ abs(z)) .* (abs(z) < pot.delta); % bad at t=0
	dpot = @(pot, z) sign(z) .* (abs(z) < pot.delta);
	shrink = @(pot, b, reg) ir_tav_shrink(b, reg, pot.delta);

% huber potential function
case 'huber'
	potk = @(pot, z) huber_pot(z, pot.delta);
	wpot = @(pot, z) huber_wpot(z, pot.delta);
	dpot = @(pot, z) huber_dpot(z, pot.delta);
	shrink = @(pot, b, reg) ir_huber_shrink(b, reg, pot.delta);

% cauchy penalty: d^2 / 2 * log(1 + (t/d)^2) (not convex!)
case 'cauchy'
	potk = @(pot, z) pot.delta.^2 / 2 .* log(1 + abs(z ./ pot.delta).^2);
	wpot = @(pot, z) 1 ./ (1 + abs(z ./ pot.delta).^2);
	dpot = @(pot, z) z ./ (1 + abs(z ./ pot.delta).^2);
	shrink = @(pot, b, reg) cauchy_shrink(b, reg, pot.delta);

% Geman&McClure penalty: d^2 / 2 * |z/d|^2 / (1 + |z/d|^2)
% Not convex!
case 'geman&mcclure'
	potk = @(pot, z) pot.delta.^2 / 2 .* abs(z/pot.delta).^2 ./ (1 + abs(z ./ pot.delta).^2);
	wpot = @(pot, z) 1 ./ (1 + abs(z ./ pot.delta).^2).^2;
	dpot = @(pot, z) z ./ (1 + abs(z ./ pot.delta).^2).^2;

% gf1: Generalized Fair 1st-order
% wpot(z) = (1 + a * |z/d|) / (1 + b * |z/d|)
case 'gf1'
	potk = @(pot, z) gf1_potk(z, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, z) (1 + pot.param(1) .* abs(z ./ pot.delta)) ...
		./ (1 + pot.param(2) .* abs(z ./ pot.delta));
	shrink = @(pot, b, reg) ...
		ir_gf1_shrink(b, reg, pot.delta, pot.param(1), pot.param(2));


% hyperbola penalty: d^2 * [ sqrt(1 + (z/d)^2) - 1 ]
case 'hyper2'
	potk = @(pot, z) pot.delta.^2 .* (sqrt(1 + abs(z ./ pot.delta).^2) - 1);
	wpot = @(pot, z) 1 ./ sqrt(1 + abs(z ./ pot.delta).^2);
	dpot = @(pot, z) z ./ sqrt(1 + abs(z ./ pot.delta).^2);

case 'hyper'
	error 'use "cauchy" or "hyper3" not "hyper" now'

% Lange1 penalty
case 'lange1'
	potk = @(pot, z) abs(z).^2 / 2 ./ (1+abs(z./pot.delta));
	wpot = @(pot, z) (1 + abs(z ./ pot.delta) / 2) ./ (1 + abs(z ./ pot.delta)).^2;

% Lange3 penalty
case {'lange3', 'fair'}
	potk = @(pot, z) pot.delta.^2 .* (abs(z./pot.delta) - log(1+abs(z./pot.delta)));
	wpot = @(pot, z) 1 ./ (1 + abs(z ./ pot.delta));
	dpot = @(pot, z) z ./ (1 + abs(z ./ pot.delta));
	% caution: no built-in shink, use 'fair-l1' if you want shrink!

% Fair potential "rounded corner" approximation to l1
case 'fair-l1'
	potk = @(pot, z) abs(z) - pot.delta .* log(1+abs(z./pot.delta));
	wpot = @(pot, z) 1 ./ (pot.delta + abs(z));
	dpot = @(pot, z) z ./ (pot.delta + abs(z));
	shrink = @(pot, b, reg) fair_l1_shrink(b, reg, pot.delta);

% li98cfs
case 'li98cfs'
	% f = @(x) atan(x) / x - 0.5; fsolve(f, 2.3)
	delta = delta / 2.3311;
	potk = @(pot, z) ir_li98cfs_potk(z, pot.delta);
	wpot = @(pot, z) ir_li98cfs_wpot(z, pot.delta);

% qgg2: q-generalized gaussian for p=2, due to Thibault, Sauer, Bouman
% q = "param", same as lange1 when q=1
case 'qgg2'
	potk = @(pot, z) z.^2 / 2 ./ (1+abs(z./pot.delta).^(2-pot.param));
	wpot = @(pot, z) (1 + abs(z ./ pot.delta).^(2-pot.param) * pot.param ...
		 / 2) ./ (1 + abs(z ./ pot.delta).^(2-pot.param)).^2;

% genhub : generalized Huber (switch between two generalized gaussians)
% same as Huber when p=2 and q=1
% p is power near 0, q is asymptotic power
case 'genhub'
	potk = @(pot, z) ir_genhub_potk(z, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, z) ir_genhub_wpot(z, pot.delta, pot.param(1), pot.param(2));

% from stevenson:94:dpr
% p = param(1), q = param(2), same as Huber when p=2 and q=1 ???
case 'stevenson94dpr'
	warn('Potential %s does not work well; use "genhub" instead', ptype)
	potk = @(pot, z) ir_stevenson94dpr_potk(z, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, z) ones(size(z), class(z)); % bogus attempt at upper bound

% tabulate derivative of pot at z_k = dz * k for k=1,...,K
% param is samples of derivative dpot(z_k)
% dpot(0) = 0, and use sample-and-hold interpolation of dpot()
% except use linear interpolation over [0 dz]
case 'table0'
	if ~isnumeric(param) || numel(param) < 10
%	if ~iscell(param) || numel(param) ~= 2
%		fail 'for table0 param must be {dz, dpot([1:K]*dz}'
		fail 'for table0 param must be [dz, dpot([1:K]*dz)]'
	end
%	param = ir_table0_setup(param{1}, col(param{2}));
	param = ir_table0_setup(param(1), col(param(2:end)));
	potk = @(pot, z) ir_table0_potk(z, pot.param);
	wpot = @(pot, z) ir_table0_wpot(z, pot.param);
	shrink = @(pot, b, reg) ir_table0_shrink(b, reg, pot.param);

% tabulate derivative of pot at z_k = dz * k for k=1,...,K
% param is samples of derivative dpot(z_k)
% dpot(0 = 0, and use linear interpolation of dpot(t)
case 'table1'
	if ~isnumeric(param) || numel(param) < 10
%	if ~iscell(param) || numel(param) ~= 2
%		fail 'for table1 param must be {dz, dpot([1:K]*dz}'
		fail 'for table1 param must be [dz, dpot([1:K]*dz)]'
	end
%	param = ir_table1_setup(param{1}, col(param{2}));
	param = ir_table1_setup(param(1), col(param(2:end)));
	potk = @(pot, z) ir_table1_potk(z, pot.param);
	wpot = @(pot, z) ir_table1_wpot(z, pot.param);
	shrink = @(pot, b, reg) ir_table1_shrink(b, reg, pot.param);

otherwise
	fail('Unknown potential "%s"', ptype)
end

if isempty(dpot) % default is z * wpot(z)
	dpot = @(pot, z) z .* wpot(pot, z);
end
if isempty(shrink)
	shrink = @ir_potential_fun_shrink;
end

end % ir_potential_fun_parse()


% ir_potential_fun_shrink()
% find argmin_z 1/2 |z - b|^2 + reg * pot(z)
% default method uses fzero() which will be slow!
function out = ir_potential_fun_shrink(pot, b, reg)
out = zeros(size(b), class(b));
opt = optimset('tolx', 1e-7);
for ii=1:numel(b)
	a = abs(b(ii));
	s = sign(b(ii));
	cost = @(z) 0.5 * (z - a).^2 + reg * pot.potk(z);
%	dfun = @(z) z - a + reg * pot.dpot(z);
	try
	%	out(ii) = s * fzero(dfun, a);
	%	out(ii) = s * fminsearch(cost, a, opt);
		z0 = s * fminsearch(cost, 0, opt);
		z1 = s * fminsearch(cost, a, opt);
		if cost(z0) < cost(z1)
			out(ii) = z0;
		else
			out(ii) = z1;
		end
	catch
		whos
		keyboard
	end
end % for
end % ir_potential_fun_shrink


% ir_table0_setup()
% dz	[1]	z spacing
% dk	[K]	dpot([1:K] * dz)
function out = ir_table0_setup(dz, dk)
sk = dk(1)*dz/2 + dz * [0; cumsum(dk(1:end-1))]; % [K]
out.dz = dz;
out.dk = dk; % [1:K] samples of dpot
out.sk = sk; % [K] cumulative sums for pot(z)
K = numel(dk);
out.K = K;
end % ir_table0_setup


% ir_table0_potk()
function out = ir_table0_potk(z, param)
dz = param.dz;
dk = param.dk;
sk = param.sk;
z = abs(z);
k = floor(z / dz);
k = min(k, param.K); % at most K
big = z > dz;
out = zeros(size(z), class(z));
out(~big) = 0.5 * dk(1) / dz * z(~big).^2;
k = k(big);
out(big) = sk(k) + dk(k) .* (z(big) - k * dz);
end % ir_table0_potk


% ir_table0_wpot()
function out = ir_table0_wpot(z, param)
dz = param.dz;
dk = param.dk; % [K]
z = abs(z);
k = floor(z / dz);
k = min(k, param.K); % at most K
out = repmat(dk(1) / dz, size(z));
big = k > 0;
k = k(big);
z = z(big);
out(big) = dk(k) ./ z;
end % ir_table0_wpot


% ir_table0_shrink()
function out = ir_table0_shrink(b, reg, param)
K = param.K;
dk = param.dk;
dz = param.dz;
zk = (1:K)' * dz;
bk = zk + reg * dk; % must be done here because depends on reg
ck = bk + dz;
tmp = [0; col([bk ck]')];
out = [zk zk+dz];
out = [0; col([zk zk+dz]')];
out = sign(b) .* interp1(tmp, out, abs(b), 'linear', 'extrap');
end % ir_table0_shrink


% ir_table1_setup()
% dz	[1]	z spacing
% dk	[K]	dpot([1:K] * dz)
function out = ir_table1_setup(dz, dk)
K = numel(dk);
k = [0:K]'; % [K+1]
zk = dz * k; % [K+1]
dk0 = [0; dk]; % [K+1] prepend sample at zero
ck = [diff(dk0); 0] / dz; % [K+1] curvatures 0:K
tmp = (dk0 - zk .* ck) * dz ...
	+ dz^2/2 * ck .* ((k+1).^2 - k.^2); % [K+1]
sk = cumsum( [0; tmp] );
out.dz = dz;
out.dk = dk; % [1:K] samples of dpot
out.ck = ck; % [K+1] curvatures 0:K
out.sk = sk; % [K] cumulative sums for pot(t)
out.K = K;
end % ir_table1_setup


% ir_table1_potk()
function out = ir_table1_potk(z, param)
sk = param.sk;
dz = param.dz;
ck = param.ck;
dk0 = [0; param.dk]; % [K+1] prepend sample at zero
z = abs(z);
k = floor(z / dz);
k = min(k, param.K); % at most K
sk = sk(1+k); % matlab indexing
ck = ck(1+k); % matlab indexing
dk = dk0(1+k);
zk = k * dz;
ck = reshape(ck, size(z));
dk = reshape(dk, size(z));
sk = reshape(sk, size(z));
out = sk + (dk - zk .* ck) .* (z - zk) + ck / 2 .* (z.^2 - zk.^2);
end % ir_table1_potk


% ir_table1_wpot()
function out = ir_table1_wpot(z, param)
dz = param.dz;
ck = param.ck;
dk = param.dk;
z = abs(z);
k = floor(z / dz);
k = min(k, param.K); % at most K
out = repmat(dk(1) / dz, size(z));
big = k > 0;
k = k(big);
z = z(big);
dk = reshape(dk(k), size(z));
ck = reshape(ck(k+1), size(z));
out(big) = (dk + (z - k * dz) .* ck) ./ z;
end % ir_table1_wpot


% ir_table1_shrink()
% unfortunately this routine works only for a single scalar reg value.
function out = ir_table1_shrink(z, reg, param)
K = param.K;
zk = (1:K)' * param.dz;
bk = zk + reg * param.dk; % must be done here because depends on reg
out = sign(z) .* interp1([0; bk], [0; zk], abs(z), 'linear', 'extrap');
end % ir_table1_shrink


% gf1_potk()
% gf1: generalized fair 1st-order potential
function pot = gf1_potk(z, delta, a, b)
atd = abs(z ./ delta);
pot = delta.^2 ./ (2 * b.^3) * ...
	(2 * b.^2 .* atd + a .* b.^2 .* atd.^2 ...
	- 2 * a .* b .* atd + 2 * (a-b) .* log(1 + b .* atd));

if 0 % symbolic check
	syms x
	syms a positive
	syms b positive
	syms z positive
	int(x*(1+a*x) / (1+b*x), x, 0, z)
end
end % gf1_potk



% ir_broken_shrink()
function out = ir_broken_shrink(z, reg, delta)
out = z ./ (1 + reg);
big = delta * (1 + reg) < abs(z);
out(big) = z(big);
end % ir_broken_shrink


% cauchy_shrink()
function out = cauchy_shrink(z, reg, delta)
z = z(:);
coef = ones(numel(z),1);
coef = [1/delta^2*coef -z/delta^2 (1+reg)*coef -z];
out = zeros(size(z));
for ii=1:numel(z)
	tmp = roots(coef(ii,:));
	pick = tmp == real(tmp); % empirically, 3rd root is often real
	if sum(pick) ~= 1 % if multiple real roots, empirically pick largest
		pick = imax(abs(tmp));
	end
	out(ii) = tmp(pick);
end
end % cauchy_shrink


% ir_huber_shrink()
function out = ir_huber_shrink(z, reg, delta)
out = z ./ (1 + reg);
big = delta .* (1 + reg) < abs(z);
if numel(reg) > 1, reg = reg(big); end
if numel(delta) > 1, delta = delta(big); end
z = z(big);
out(big) = z .* (1 - reg .* delta ./ abs(z));
end % ir_huber_shrink


% fair_l1_shrink()
function out = fair_l1_shrink(z, reg, delta)
out = sign(z) .* (abs(z) - (delta + reg) ...
	+ sqrt( (delta + reg - abs(z)).^2 + 4 * delta .* abs(z) )) ./ 2;
end % fair_l1_shrink


% ir_gf1_shrink()
function out = ir_gf1_shrink(z, reg, delta, a, b)
u = a / delta;
v = b / delta;
out = sign(z) .* (v .* abs(z) - (1+reg) ...
	+ sqrt( (1 + reg - v.*abs(z)).^2 + 4 * (v + reg .* u) .* abs(z) )) ...
	./ (2 * (v + reg .* u));
end % ir_gf1_shrink


% ir_lpp_1_2_shrink()
% for l_p with p=1/2
% for a > 0 and z = t^2 the minimizer solves t^3 - a t + reg/2 = 0
function out = ir_lpp_1_2_shrink(b, reg)
sb = sign(b); a = abs(b);
%reg = 1;
%a = 50; roots([1 0 -a reg/2])
%a = linspace(0, 4, 401);
out = zeros(size(b), class(b));
%big = true(size(b));
%{
% loop method
for ib=1:numel(b)
	tmp = roots([1 0 -a(iz) reg/2])'
	out(ib) = max(tmp).^2;
end
%}
%{
% hyperbolic method for one real root:
big = 27 * (reg/2)^2 > 4 * a.^3;
t0 = -2 * sqrt(a/3) .* cosh(1/3 * acosh(3/2 * (reg/2) ./ a .* sqrt(3./a)));
%}
%{
% https://en.wikipedia.org/wiki/Cubic_function#Vieta.27s_substitution
% vieta solves t^3 + p t + q = 0 using t = w - p / (3w), i.e. p=-a and q=reg/2
% for which w^6 + q w^3 - p^3/27 = 0
% i.e. (w^3)^2 + reg/2 (w^3) + a^3/27 = 0
w3 = (-reg/2 + sqrt((reg/2).^2 - 4 * a.^3/27)) / 2; % w^3 solution by quad form
w1 = w3 .^ (1/3);
t0 = w1 + a ./ (3*w1);
out(big) = t0(big).^2;
%}
% trigonometric method for three real roots:
% https://en.wikipedia.org/wiki/Cubic_function#Trigonometric_method_for_three_real_roots
t0 = 2 * sqrt(a/3) .* cos(1/3 * acos(3*(reg/2)/2./(-a) .* sqrt(3./a)));
a_1_2 = 3/2 * reg^(2/3);
big = a > a_1_2;
out(big) = t0(big).^2;
% plot(z, out, '-o', z, z, ':', z, z-reg, '--')
out = out .* sb;
end % ir_lpp_1_2_shrink


% ir_lpp_4_3_shrink()
% for l_p with p=4/3
% from eqn. (4.5) in chaux:07:avf
function out = ir_lpp_4_3_shrink(b, reg)
sb = sign(b); a = abs(b);
x = sqrt(a.^2 + 256 * reg^3 / 729);
out = sb .* (a + 4 * reg / (3 * 2^1/3) * ((x - a).^1/3 - (x + a).^1/3));
end % ir_lpp_4_3_shrink


% ir_tav_shrink()
function out = ir_tav_shrink(z, reg, delta)
out = zeros(size(z), class(z));
big = reg < abs(z) & abs(z) < reg + delta;
out(big) = z(big) .* (1 - reg ./ abs(z(big)));
big = reg + delta <= abs(z);
out(big) = z(big);
end % ir_tav_shrink


% ir_li98cfs_potk()
function pot = ir_li98cfs_potk(z, d)
pot = d.^2 .* ((z ./ d) .* atan(z ./ d) - 0.5 * log(1 + (z ./ d).^2));
end % ir_li98cfs_potk


% ir_genhub_potk()
function pot = ir_genhub_potk(z, d, p, q)
pot = 0.5 * abs(z) .^ p .* (abs(z) <= d) ...
 + 0.5 * (p ./ q .* d .^ (p-q) .* abs(z) .^ q ...
 + (1 - p ./ q) .* d .^ p) .* (abs(z) > d);
end % ir_genhub_potk

% ir_genhub_wpot()
function pot = ir_genhub_wpot(z, d, p, q)
pot = p / 2 .* (d .^ (p-q)) .* (abs(z) .^ (q-2)) .* (abs(z) > d);
ii = abs(z) <= d;
pot(ii) = p / 2 .* (abs(z(ii)) .^ (p-2));
%pot = p / 2 .* abs(t) .^ (p-2) .* (abs(t) <= d) ...
% + p / 2 .* d .^ (p-q) .* abs(t) .^ (q-2) .* (abs(t) > d);
end % ir_genhub_wpot


% ir_stevenson94dpr_potk()
function pot = ir_stevenson94dpr_potk(z, d, p, q)
% pr [d p q]
tmp1 = (0.5 * abs(z) .^ p) .* (abs(z) <= d);
tmp2 = 0.5 * ( (p .* (d .^ (p-1)) .* abs(z) - p .* (d .^ p) ...
 + (1 ./ q) .^ (1 ./ (q-1)) ) .^ q ...
 + d .^ p - (1 ./ q) .^ (q ./ (q-1)) ) .* (abs(z) > d);
pot = tmp1 + tmp2;
end % ir_stevenson94dpr_potk


% ir_tonumeric()
% convert x to type of y
function out = ir_tonumeric(x, y)
switch class(y)
case 'double'
	out = double(x);
case 'single'
	out = single(x);
otherwise
	fail('unknown type %s', class(y))
end
end % ir_tonumeric


% ir_potential_fun_plot()
function dummy = ir_potential_fun_plot(pot)
z = linspace(-1,1,101)*2;
if isvar('pot.delta')
	z = z * pot.delta;
end
if ~im, return, end
im plc 2 2
im subplot 1
plot(z, pot.potk(z), '-', z, z.^2/2, ':', z, abs(z), ':')
axis([min(z) max(z) minmax(pot.potk(z))'])

im subplot 2
plot(z, pot.dpot(z), '-', z, z, ':')
axis([min(z) max(z) minmax(pot.dpot(z))'])

im subplot 3
plot(z, pot.wpot(z), '-', z, 1+0*z, ':')
axis([min(z) max(z) 0 1])

im subplot 4
reg = 1;
tmp = pot.shrink(z, reg);
plot(z, tmp, '-', z, z, ':')
axis([min(z) max(z) minmax(tmp)'])

dummy = [];
end % ir_potential_fun_plot


% ir_potential_fun_test()
% test routine
% examine potential functions after rescaling.
function ir_potential_fun_test

delta = 10; zmax = 4 * delta; reg = 1.5 * delta;
plist = potential_fun('list');
%plist = {'l1', 'fair-l1'}; delta = 0.5;
%plist = {'quad', 'li98cfs', 'hyper3', 'huber2'}; % show li98cfs roughly hyper3
%plist = {'quad', 'genhub', 'huber', 'stevenson94dpr'};
%plist = {'genhub'}
%plist = {'hyper3', 'qgg2', 'huber'}; delta = 10; zmax = 50;
%plist = {'lange3', 'qgg2'}; zmax = 200;
%plist = {'qgg2', 'gf1-fit'};
%plist = potential_fun('list1'); delta = 0.5;
%plist = {plist{:}, 'quad', 'gf1', 'huber', 'broken'}; % todo: more!
%plist = {'qgg2', 'table1', 'table0'};
%plist = {'qgg2', 'table1'}; fname = 'fig_reg_pot_table1_qgg2';
%plist = {'qgg2', 'table0'}; fname = 'fig_reg_pot_table0_qgg2';
%plist = {'qgg2', 'table1', 'gf1'}; fname = 'fig_reg_pot_table0_gf1_qgg2';
%plist = {'cauchy', 'l1'};
zz = zmax * linspace(-1, 1, 2001)';
bb = 3 * max(delta + reg, delta * (1+reg)) * linspace(-1, 1, 301)';
ps = [];
lshrink = {};
for ii=1:numel(plist)
	ptype = plist{ii}; % pr ptype
	if streq(ptype, 'quad')
		leg{ii} = ptype;
	else
		leg{ii} = [ptype sprintf(', $\\delta = %g$', delta)];
	end

	switch ptype
	case {'gf1', 'gf1-fit'}
		param = potential_fun('gf1-fit', nan, {'qgg2', [1 10], 1.2});
		ptype = 'gf1'; % trick
		leg{ii} = [leg{ii} sprintf(' %.3g %.4f', param(1), param(2))];
	case 'qgg2'
		param = 1.2;
		leg{ii} = [leg{ii} sprintf(', $q = %g$', param)];
	case 'genhub'
		param = [2.0 1.2];
		leg{ii} = [leg{ii} sprintf(', $p=%g, q=%g$', param(1), param(2))];
	case 'stevenson94dpr'
		param = [2 1.2];
		leg{ii} = [leg{ii} sprintf(', $p=%g, q=%g$', param(1), param(2))];
	case {'table0', 'table1'}
		dz = delta / 20;
		tmp = potential_fun('qgg2', delta, 1.2);
%		param = {dz, tmp.dpot([1:1e4] * dz)};
		param = [dz, tmp.dpot([1:1e4] * dz)];
%		tmp = sprintf(' K = %d dz = %g', numel(param{2}), param{1});
		tmp = sprintf(', $K = %d$, $\\Delta t = %g$', numel(param)-1, param(1));
		leg{ii} = [leg{ii} tmp];
	otherwise
		param = [];
	end

	pot = potential_fun(ptype, delta, param);
	pp(:,ii) = pot.potk(zz);
	pw(:,ii) = pot.wpot(zz);
	pd(:,ii) = pot.dpot(zz);

	% replace 0 with 1 to make (slow) figure showing fzero-based shrinkage
	if 0 || ~streq(func2str(pot.meth.shrink), 'ir_potential_fun_shrink')
		tmp = pot.shrink(bb, reg);
		if any(tmp ~= 0)
			ps(:,end+1) = pot.shrink(bb, reg);
			lshrink{end+1} = leg{ii};
		end
	end

	if 0 % test vs old
		try
			opot = potential_func(ptype, delta, param);
			opp = opot.potk(opot, zz);
			opw = opot.wpot(opot, zz);
			opd = opot.dpot(opot, zz);
		catch
			printm('%s no old version', ptype)
		%	keyboard % problem with private scope of ir_li98cfs_wpot
			continue
		end
		if ~isequal(opp, pp(:,ii)), 'p bug', ptype, keyboard, end
		if ~isequal(opw, pw(:,ii)), 'w bug', ptype, keyboard, end
		if ~isequal(opd, pd(:,ii)), 'd bug', ptype, keyboard, end
	end
end

if im
%	set(0,'DefaultAxesLineStyleOrder', '-|:')
	clf, pl = @(i) subplot(410 + i);
	pl(1), plot(zz, pp), title 'potk'
	axis tight, axisy([-0.0 2.5] * delta^2)
	ir_legend(leg, 'location', 'north')
	pl(2), plot(zz, pw), title 'wpot'
	axis tight, axisy(0, 1.1), ytick([0 1])
	pl(3), plot(zz, pd), title 'dpot'
	axis tight, axisy([-1 1] * 1.3 * delta)
	xlabelf '$z$'
%	ir_savefig eps_c fname
	% check derivatives
	pl(4)
	plot(zz, pd)
	title 'dpot check'
	hold on
	d = diffc(pp) / (zz(2)-zz(1));
	plot(zz(1:end-1), d(1:end-1,:), '--', ...
		zz(1:end-1), d(1:end-1,:)-pd(1:end-1,:), ':')
	hold off
	axis tight, axisy([-1 1] * 1.3 * delta)
end

if 1 && im
	if ~isempty(lshrink)
		prompt
		clf
		plot(bb, ps, '-', bb, bb, '-')
		axis equal, axis square
		axis([-1 1 -1 1]*400)
		xtick([-1 0 1] * 400)
		ytick([-1 0 1] * 400)
		ir_legend(lshrink, 'location', 'southeast')
		xlabelf '$c$'
%		ylabelf 'xhat(z)'
		ylabelf '$\hat{z}(c)$'
		grid
%		ir_savefig cw fig_reg_pot_table0_qgg2_shrink

		if 0 % figure for book
			plot(bb, ps(:,[3 2]) - repmat(ps(:,1), [1 2]), 'o')
			plot(	bb, ps(:,3) - ps(:,1), 'bo', ...
				bb, ps(:,2) - ps(:,1), 'gx')
			axis([0 500 -0.05 0.5]), ytick([0 0.05 0.1 0.5])
			xtick([0 500])
			ir_fontsize label 18
			ir_fontsize text 18
			xlabelf 'c', ylabelf 'shrinkage error for QGG2'
			ir_legend(lshrink([3 2]), 1)
%			ir_savefig cw fig_reg_pot_table01_qgg2_shrinker
			keyboard
		end
	end
end

% check 'scale' option
pot1 = potential_fun('lange3', 10, []);
pot2 = potential_fun('lange3', 10, [], 'scale', 2);
jf_equal(2 * pot1.potk(zz), pot2.potk(zz))
jf_equal(2 * pot1.dpot(zz), pot2.dpot(zz))
jf_equal(2 * pot1.wpot(zz), pot2.wpot(zz))
 
end % ir_potential_fun_test


% ir_potential_fun_test_lpp()
function ir_potential_fun_test_lpp
params = [0 0.01 0.1 1/2 3/4 1 4/3];
params = [0 1/2 1];
reg = 4*2^4; zmax = 2 * reg;
zz = zmax * linspace(-1, 1, 2001)';
bb = 2 * reg * linspace(-1, 1, 201);
%a_1_2 = 3 * (reg/4)^(2/3) % saddle point for p=1/2
a_1_2 = 3/2 * reg^(2/3); % breakpoint for p=1/2
if 0
	a = a_1_2;
	pot = potential_fun('lpp', nan, 1/2);
	z0 = ir_potential_fun_shrink(pot, a+1e-9, reg)
	cost = @(z) 1/2 * (z - a).^2 + reg * pot.potk(z);
	t0 = sqrt(z0)
	t0^3 - 2*a*t0 + 2*reg % should be 0
	t0^3 - a*t0 + reg/2 % should be 0
	3/2*reg - a * t0 % should be 0
	(3^(3/2) * reg)/(4*a^(3/2)) % should be 1/sqrt(2)
	cost(0)
	cost(z0)
%	syms x; solve(3 * acos(x) - acos(-x), x) % 1/sqrt(2)
return
end
tmp = [sqrt(2*reg) reg a_1_2];
tmp = outer_sum(tmp, [-1 +1]*1e-5);
bb = sort([bb col(tmp)' -col(tmp)'])'; % trick
ps = [];
lshrink = {};
if 0
	pot = potential_fun('lpp', nan, 1/2);
	potk = @(z) pot.potk(z);

%	y = sqrt(2*reg); % l0
%	y = reg; % l1
	y = a_1_2 + 1e-9; % p=1/2
	if 0
		ir_potential_fun_shrink(pot, y, reg)
		ir_potential_fun_shrink(pot, y, reg)
	return
	end

	if 0
		z1 = pot.shrink(y, reg)
		d1 = z1 - y + reg * pot.dpot(z1) % 1st derivative
		d2 = 1 + reg * (1/2)*(-1/2) * z1^(-3/2) % 2nd derivative > 0?

		dfit = 1/2 * (y - zz).^2;
		cost = dfit + reg * potk(zz);
		plot(zz, dfit, '--', zz, reg*potk(zz), '-', zz, cost, '-')
		tmp = minmax(zz(cost <= 3*min(cost)));
		tmp(1) = min(tmp(1), -10);
		tmp(2) = max(tmp(2), reg/3);
		xlim(tmp), ylim([0.0 2*1.2] * min(cost)), grid
	%	ytick([0 round(min(cost))])
		ytick([0 min(cost)])
		xtick([0 z1 a_1_2]), xlabelf '$z$', ylabelf '$\Psi(z)$'
	return
	end

	zs = ir_potential_fun_shrink(pot, bb, reg);
	zp = pot.shrink(bb, reg);

	plot(bb, zs, '-o', bb, zp, '.-', bb, bb, ':')
	legend({'fmin', 'trig'}, 'location', 'northwest')
	xlabelf '$b$'
	ylabelf '$\hat{z}(b)$'
	axis equal, axis square
	axis([-1 1 -1 1]*1.1*reg)
	tmp = [sqrt(2*reg) a_1_2 reg]; % l0 and l1 breakpoints
	xtick([-tmp 0 tmp])
	ytick([-1 0 1] * reg)
	grid
%	ir_savefig cw fig_?
return
end

for ip=1:length(params)
	pot = potential_fun('lpp', 0, params(ip));

	% replace 0 with 1 to make (slow) figure showing fzero-based shrinkage
	if 0 || ~streq(func2str(pot.meth.shrink), 'ir_potential_fun_shrink')
		tmp = pot.shrink(bb, reg);
		if any(tmp ~= 0)
			ps(:,end+1) = pot.shrink(bb, reg);
			lshrink{end+1} = sprintf('p=%g', params(ip));
		end
	end

	clf
	plot(bb, ps, '-', bb, bb, ':')
	legend(lshrink{:}, 'location', 'southeast')
	xlabelf '$b$'
	ylabelf '$\hat{z}(b)$'
	axis equal, axis square
	axis([-1 1 -1 1]*1.5*reg)
	tmp = [sqrt(2*reg) a_1_2 reg]; % l0 and l1 breakpoints
	xtick([-tmp 0 tmp])
	ytick([-1 0 1] * reg)
	grid
%	ir_savefig cw fig_?
end

end % ir_potential_fun_test_lpp()
