 function proj = cylinder_proj(cg, params, varargin)
%function proj = cylinder_proj(cg, params, [options])
%|
%| Compute set of 2d line-integral projection views of (elliptical) cylinder(s).
%| Works for these 3D geometries:
%|	parallel beam
%|	flat-detector cone-beam
%|	arc-detector cone-beam (3rd generation CT)
%|	todo: would be nice to have tent cone-beam too!
%|
%| in
%|	cg			ct_geom()
%|	params [ne 8]		elliptical cylinder parameters:
%|		[centx centy centz  radx rady zlength  angle_degrees  amplitude]
%|
%| options
%|	oversample		oversampling factor for emulating "strips"
%|				(to account for finite detector size)
%|
%| out
%|	proj	[ns nt na]	projection views
%|
%| Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| Rebecca Malinas, Ajay Paidi and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test'), ir_cylinder_proj_test, return, end
if nargin < 2, ir_usage, end

warn 'todo: currently incorrect for certain angles!'

% defaults
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = ir_cylinder_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % cylinder_proj()


% ir_cylinder_proj_do()
function proj = ir_cylinder_proj_do(params, ss, tt, ...
		beta, ... % [na] [radians]
		source_zs, dso, dod, dfs, oversample)

if size(params, 2) ~= 8, fail '8 parameters per cylinder', end 

if oversample > 1
	ds = ss(2) - ss(1);
	dt = tt(2) - tt(1);
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6)
		fail 'uniform spacing required for oversampling'
	end
	No = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No-1):2:(No-1)]'/(2*No)*ds, ss(:)'); % [No ns]
	tt = outer_sum([-(No-1):2:(No-1)]'/(2*No)*dt, tt(:)'); % [No nt]
	proj = ir_cylinder_proj_work(params, ss(:), tt(:), beta, source_zs, ...
			dso, dod, dfs);
	proj = downsample3(proj, [No No 1]);
else
	proj = ir_cylinder_proj_work(params, ss(:), tt(:), beta, source_zs, ...
			dso, dod, dfs);
end

end % ir_cylinder_proj_do()


% ir_cylinder_proj_work()
function proj = ir_cylinder_proj_work(params, ss, tt, ...
		beta, ... % [radians]
		source_zs, dso, dod, dfs, oversample)

% determine equivalent parallel-beam projection coordinates, at beta=0
ns = numel(ss);
nt = numel(tt);
[sss, ttt] = ndgrid(ss, tt);

if isinf(dso) % parallel beam
	uu = sss; % transaxial
	vv = ttt; % axial
	azim0 = zeros(size(uu)); % azimuthal for starting angle
	polar = zeros(size(uu));

elseif isinf(dfs) % cone-beam with flat detector
	[uu, vv, azim0, polar] = ir_coord_cb_flat_to_par(sss, ttt, dso, dod);

elseif dfs == 0 % cone-beam with (3rd-gen) arc detector
	[uu, vv, azim0, polar] = ir_coord_cb_arc_to_par(sss, ttt, dso, dod);

else
	fail 'not done'
end

clear sss ttt % using parallel-beam coord. hereafter

cpolar = cos(polar);
spolar = sin(polar);
proj = zeros(ns, nt, numel(beta), 'single');

% loop over cylinders
for ip = 1:size(params,1)
	par = params(ip,:);

	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	zh = par(6) / 2; % half of z_length of cylinder
	eang = deg2rad(par(7)); % xy-plane rotation
	val = par(8);

	for ib = 1:numel(beta)
		az = beta(ib) + azim0;

		% shift property of 3D transform:
		cz_eff = cz - source_zs(ib); % center relative to source
		ushift = cx * cos(az) + cy * sin(az);
		vshift = (cx * sin(az) - cy * cos(az)) .* spolar + cz_eff * cpolar;
		az = az - eang;

		% substitute parametric equation of ray (p = p0 + l*p1)
		% in equation for infinite elliptical cylinder

		p1 = (uu-ushift) .* cos(az) + (vv-vshift) .* sin(az) .* spolar;
		p2 = (uu-ushift) .* sin(az) - (vv-vshift) .* cos(az) .* spolar;
		p3 = (vv-vshift) .* cpolar;

		e1 = -sin(az) .* cpolar; % direction cosines of ray
		e2 = cos(az) .* cpolar;
		e3 = spolar;

		A = e1.^2 / rx^2 + e2.^2 / ry^2;
		B = p1 .* e1 / rx^2 + p2 .* e2 / ry^2;
		C = p1.^2 / rx^2 + p2.^2 / ry^2 - 1;

		% calculate l at intersection points of ray with inf. cylinder
		det = B.^2 - A .* C;
		good = det >= 0; % real roots => ray intersects inf. cylinder
		tmp = sqrt(det(good));
		l0 = zeros(size(det), 'single');
		l1 = zeros(size(det), 'single');
		A = A(good);
		B = B(good);
		l0(good) = (-B - tmp) ./ A;
		l1(good) = (-B + tmp) ./ A;

		% z values at the points of intersection
		z0 = p3 + l0 .* e3;
		z1 = p3 + l1 .* e3;

		% re-arrange so that z0 <= z1
		zswap = z0 > z1;
		[z0, z1] = ir_swap(z0, z1, zswap);

		% truncate to (usually finite) extent of cylinder
		zmin = max(z0, -zh);
		zmax = min(z1, zh);

		% scale line-integral by fraction of z-range within cylinder
		l_int = zeros(size(l1), 'single');
		zok = good & (z1 ~= z0);
%		l_int = abs(l1 - l0) .* max(zmax - zmin, 0) ./ abs(z1 - z0);
		tmp0 = abs(l1 - l0);
		tmp1 = tmp0 .* max(zmax - zmin, 0);
		tmp2 = abs(z1 - z0);
		l_int(zok) = tmp1(zok) ./ tmp2(zok);

		zeq = good & (z1 == z0) & (-zh < z0) & (z0 < zh); % e3=0 rays
		l_int(zeq) = tmp0(zeq);

		proj(:,:,ib) = proj(:,:,ib) + val * l_int;
	end % ib
end % ip

%{ old way

%		% re-arrange l0 and l1 according to new z ordering
%		[l0 l1] = ir_swap(l0, l1, zswap);
%
%		if any(col(l0 > l1)), fail('l0 > l1'), end

		t0 = zeros(size(l0), 'single');
		t1 = zeros(size(l1), 'single');

		% store l-values for cases entirely within finite cylinder
		tmp = (-zh < z0) & (z0 < zh) & (-zh < z1) & (z1 < zh);
		t0(tmp) = l0(tmp);
		t1(tmp) = l1(tmp);

		% set upper and lower z-bounds on the infinite z cylinder
		l_up = (zh - p3) ./ e3; % bug: divide by 0
		l_down = (-zh - p3) ./ e3;

		% if z values exceed this bound, then set those z values to the
		% limiting value
		tmp = (zh < z1) & (-zh < z0) & (z0 < zh); % "index_up"
		t0(tmp) = l0(tmp);
		t1(tmp) = l_up(tmp);

		tmp = (z0 < -zh) & (-zh < z1) & (z1 < zh); % "index_down"
		t0(tmp) = l_down(tmp);
		t1(tmp) = l1(tmp);

 		% ray entering and leaving cylinder at opposite ends:
		tmp = (z0 < -zh) & (zh < z1); % "index_both"
		t0(tmp) = l_down(tmp);
		t1(tmp) = l_up(tmp);

		proj(:,:,ib) = proj(:,:,ib) + val * abs(t1 - t0);
%}

end % ir_cylinder_proj_work()


% ir_swap()
% swap specified elements
function [o0, o1] = ir_swap(i0, i1, swap)
o0 = i0;
o1 = i1;
o0(swap) = i1(swap);
o1(swap) = i0(swap);

end % ir_swap()


% ir_cylinder_proj_test()
% internal test routine
function ir_cylinder_proj_test

param = [ ... % defrise phantom type disks
%	[20 0 -82	80 80 20 0 10];
	[20 0 -42	80 80 20 0 10];
	[20 0 +0	80 80 20 0 10];
	[20 0 42	80 80 20 0 10];
%	[20 0 82	80 80 20 0 10];
];

param = [20 10 -12	180 70 80 0 10]; % quick simple cylinder test

param = [20 10 -12	200 200 800 0 10]; % quick simple cylinder test
% todo

printm('max line integral roughly %g', ...
	max(param(4:5))*2*param(8))
%	max(param(4)*abs(ig.dx), param(5)*abs(ig.dy))*2*param(8))

fun_proj = @(cg) cylinder_proj(cg, param, 'oversample', 3); % analytical
fun_im = @(ig) cylinder_im(ig, param, 'oversample', 2, 'checkfov', true);

ir_proj3_compare1(fun_proj, fun_im, 'chat', 1, 'dsd', 949, 'dfs', 0, ...
	'arg_ct_geom', {'na', 1}, ...
	'nt', 125-0, ... % stress odd sin(polar) = 0
	'downp', 2, ...
	'source_z0', 1*7, 'pitch', 0);

return

ir_proj3_compare1(fun_proj, fun_im, 'downi', 4, ...
	'chat', 1, 'dsd', inf, 'dfs', 0, 'pitch', 0);

ir_proj3_compare1(fun_proj, fun_im, 'chat', 1);


down = 4;
oversample = 1;

im plc 2 2

dfs_list = [0 inf];
for ii = 1:numel(dfs_list)
	dfs = dfs_list(ii);
	cg = ct_geom('fan', 'ns', 272, 'nt', 256, 'na', 112, ...
		'ds', 4, 'dsd', 949, 'dod', 408, 'dfs', dfs, ...
		'offset_s', 0.25, 'down', down);

	proj = cylinder_proj(cg, param, 'oversample', oversample);
	im(ii, proj, 'matlab cone-beam projections'), cbar
	titlef('matlab cone-beam projections, Dfs=%g', dfs)
	im(2+ii, proj(:,:,end/2))
end
%im(proj(:,:,24))

end % ir_cylinder_proj_test()
