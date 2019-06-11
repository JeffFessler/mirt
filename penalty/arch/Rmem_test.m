% Rmem_test.m

% check match between Reg1 and Rmem
if 1 | ~isvar('Rm'), printm 'Rm'
	if 0
		ig = image_geom('nx', 512, 'ny', 496, 'nz', 16, 'fov', 500, ...
			'down', 4);
	else % 2d
		ig = image_geom('nx', 100, 'ny', 80, 'dx', 4.2);
	end
	ig.mask = ig.circ(200) > 0;
	kappa = 2 + 1 * cos(2*pi*ig.x/ig.fov) * sin(2*pi*ig.y/ig.fov)';
	kappa = single(kappa) .* ig.mask;
%	kappa = repmat(kappa, [1 1 ig.nz]) .* ig.mask;
%	kappa = ig.mask;

	f.l2b = 3;
%	f.l2b = 0;
	f.pot = 'hyper3'; f.delta = 30; f.pot_arg = {f.pot, f.delta};
%	f.pot = 'quad'; f.delta = inf; f.pot_arg = {f.pot};
	f.offsets = '3d:26';
	f.offsets = [ig.nx+1];
	f.offsets = '';
	Rg = Reg1(kappa, 'pot_arg', f.pot_arg, 'beta', 2^f.l2b, ...
		'offsets', f.offsets, 'type_penal', 'mat');
%		'offsets', f.offsets, 'type_penal', 'mex');
%		'edge_type', 'tight', ...
	Rm = Rmem(kappa, 'pot_arg', f.pot_arg, 'beta', 2^f.l2b, ...
		'distance_power', Rg.distance_power, ...
		'offsets', f.offsets);
%		'edge_type', 'tight', ...

	x = ig.unitv;

	if 1 % check C1
		jf_equal(Rm.C1 * x, Rg.C1 * x)
		jf_equal(Rm.C1 * x(ig.mask), Rg.C1 * x(ig.mask))
	end

	if 1 % check penalty value
		jf_equal(Rg.penal(Rg, x), Rg.penal(Rg, x(ig.mask)))
		jf_equal(Rm.penal(Rm, x), Rm.penal(Rm, x(ig.mask)))
		equivs(Rg.penal(Rg, x), Rm.penal(Rm, x))
	end

	if 1 % check cgrad
		g1 = Rm.cgrad(Rm, x);
		g2 = ig.embed(Rm.cgrad(Rm, x(ig.mask)));
		jf_equal(g1, g2)
		g2 = Rg.cgrad(Rg, x);
im pl 2 2, im(1, g1), im(2, g2), im(3, g1-g2), im(4, ig.mask)
		equivs(g1, g2)
		g3 = ig.embed(Rg.cgrad(Rg, x(ig.mask)));
		jf_equal(g2, g3)
		clear g1 g2 g3
	end

	if 0 % check denom
	return
	end

	%d = Rm.diag(R);
prompt
end

% Robject designed to match Rmem'
if 1 | ~isvar('Ro'), printm 'Ro'
	tmp = repmat(kappa.^2, [1 1 1 length(Rm.offsets)]);
	Ro = Robject(ig.mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
		'offsets', Rm.offsets, ...
		'potential', f.pot, 'delta', f.delta, ...
		'type_denom', 'matlab', ...
		'distance_power', Rg.distance_power, 'user_wt', tmp);
	clear tmp
end

if 1
	cpu etic
	g1 = Rm.cgrad(Rm, x);
	cpu etoc 'Rm cgrad time'
	cpu etic
	g2 = Ro.cgrad(Ro, x(ig.mask));
	g2 = embed(g2, ig.mask);
	cpu etoc 'Ro cgrad time'
	equivs(g1, g2)
end

if 0 % these tests superceded by Reg1_test.m
	tmp = sprintf('cgrad%d,offset', Rm.order);
	beta = 2^f.l2b ...
		./ penalty_distance(Rm.offsets(:), Rm.dim) .^ Rm.distance_power;
	g31 = penalty_mex(tmp, single(x), single(kappa), ...
		int32(Rm.offsets), single(beta), ...
		Rm.pot_arg{1}{1}, single(Rm.pot_arg{1}{2}), int32(1));
	g34 = penalty_mex(tmp, single(x), single(kappa), ...
		int32(Rm.offsets), single(beta), ...
		Rm.pot_arg{1}{1}, single(Rm.pot_arg{1}{2}), int32(4));
	jf_equal(g31, g34)
	max_percent_diff(g1, g31)
	max_percent_diff(g2, g31)
end

if 1
	% compare them
	%x = ig.unitv;
	pr Rm.penal(Rm, x)
	pr Ro.penal(Ro, x(ig.mask))
	max_percent_diff 'Rm.penal(Rm,x)' 'Ro.penal(Ro,x(ig.mask))'
end

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = mask;
rng(0)
x = rand(size(mask));
x = x(mask);
max_percent_diff(R0.penal(R0, x), R1.penal(R1, x))
max_percent_diff(R0.diag(R0), R1.diag(R1))
max_percent_diff(R0.cgrad(R0, x), R1.cgrad(R1, x))
max_percent_diff(R0.denom(R0, x), R1.denom(R1, x))
if 0
	d0 = R0.denom(R0, x);
	d1 = R1.denom(R1, x);
	d0 = embed(d0, mask);
	d1 = embed(d1, mask);
	im clf
	im(221, d0)
	im(222, d1)
	im(223, d1-d0)
end
