 function proj = cuboid_proj(cg, params, varargin)
%function proj = cuboid_proj(cg, params, varargin)
%|
%| Compute set of 2d line-integral projection views of cuboids.
%| Works for both parallel-beam and cone-beam geometry.
%|
%| in
%|	cg			ct_geom()
%|	params [ne 9]		cuboid parameters:
%|			[x_center y_center z_center  x_radius y_radius z_radius
%|				xy_angle_degrees z_angle_degrees  amplitude]
%| options
%|	oversample		over-sampling factor for emulating "strips"
%|				(to account for finite detector size)
%|
%| out
%|	proj	[ns nt na]	projection views
%|
%| Yong Long, 2008-08-28, University of Michigan, adapted from ellipsoid_proj()
%| 2013-05-24 Jeff Fessler

if nargin == 1 && streq(cg, 'test'), cuboid_proj_test, return, end
if nargin < 2, ir_usage, end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = cuboid_proj_do(params, cg.s, cg.t, cg.dt, cg.ar, cg.source_zs, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % cuboid_proj()


% cuboid_proj_line1()
function [lxmin lxmax] = cuboid_proj_line1(rx, p1, e1)
tmp = (e1 == 0);
e1(tmp) = inf;
% bounds of l corresponding to rect(x/rx)
lxmin = (-rx/2 - p1) ./ e1;
lxmax = ( rx/2 - p1) ./ e1;
% re-arrange the bounds so that lxmin contains the minimum l values
% and lxmax contains the maximum l values
temp = lxmin;
lxmin = min(lxmin, lxmax);
lxmax = max(temp, lxmax);
% exclude points where e1=0 by setting lxmin = -Inf and lxmax = Inf
lxmin(tmp) = -inf;
lxmax(tmp) = inf;
end % cuboid_proj_line1()


% cuboid_proj_do()
function proj = cuboid_proj_do(params, ss, tt, dt, ...
		beta, ... % [radians]
		source_zs, dso, dod, dfs, oversample)

if size(params, 2) ~= 9, error '9 parameters per cuboid', end

if oversample > 1
	ds = ss(2) - ss(1);
%	dt = tt(2) - tt(1); % fails if nt=1 so pass dt instead
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6)
		fail 'uniform spacing required for oversampling'
	end
	No = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No-1):2:(No-1)]'/(2*No)*ds, ss(:)'); % [No ns]
	tt = outer_sum([-(No-1):2:(No-1)]'/(2*No)*dt, tt(:)'); % [No nt]
	proj = cuboid_proj_do(params, ss(:), tt(:), dt/No, ...
		beta, source_zs, dso, dod, dfs, 1);
	proj = downsample3(proj, [No No 1]);
return
end


% determine equivalent parallel-beam projection coordinates, at beta=0
ns = length(ss);
nt = length(tt);
[sss ttt] = ndgrid(ss, tt);

if isinf(dso) % parallel beam
	uu = sss;
	vv = ttt;
	azim0 = zeros(size(sss));
	polar = zeros(size(sss));

elseif isinf(dfs) % cone-beam with flat detector
	[uu vv azim0 polar] = ir_coord_cb_flat_to_par(sss, ttt, dso, dod);

elseif dfs == 0 % cone-beam with arc detector
	[uu vv azim0 polar] = ir_coord_cb_arc_to_par(sss, ttt, dso, dod);

else
	fail 'not done'
end

clear sss ttt

cpolar = cos(polar);
spolar = sin(polar);
proj = zeros(ns, nt, numel(beta));

% loop over cuboids
for ip = 1:size(params,1)
	par = params(ip,:);

	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	rz = par(6);
	xang = deg2rad(par(7)); % xy-plane rotation
	zang = deg2rad(par(8)); % z-plane rotation
	if zang, error 'z rotation not done', end
	val = par(9);

	for ib = 1:length(beta)
		az = beta(ib) + azim0;

		% shift property of 3D transform:
		cz_eff = cz - source_zs(ib); % center relative to source
		ushift = cx * cos(az) + cy * sin(az);
		vshift = (cx * sin(az) - cy * cos(az)) .* spolar + cz_eff * cpolar;

		az = az - xang;
		us = uu - ushift;
		vs = vv - vshift;
		p1 = us .* cos(az) + vs .* sin(az) .* spolar;
		p2 = us .* sin(az) - vs .* cos(az) .* spolar;
		p3 = vs .* cpolar;

		e1 = -sin(az) .* cpolar; % x = p1 + l*e1
		e2 = cos(az) .* cpolar; % y = p2 + l*e2
		e3 = spolar; % z = p3 + l*e3

		[lxmin lxmax] = cuboid_proj_line1(rx, p1, e1);
%		clear p1
		[lymin lymax] = cuboid_proj_line1(ry, p2, e2);
%		clear p2
		[lzmin lzmax] = cuboid_proj_line1(rz, p3, e3);
%		clear p3

		lmin = max(lxmin, lymin);
		lmin = max(lmin, lzmin); % lower bound for l

		lmax = min(lxmax, lymax);
		lmax = min(lmax, lzmax); % upper bound for l

		ll = max(lmax - lmin, 0); % intersection only when lmax > lmin

		% cases where e(k) = 0 (rays along axes)

		tmp = e1 == 0; % sin(az) = 0
		zero_e = (-rx/2 <= us(tmp)) & (us(tmp) <= rx/2);
		ll(tmp) = ll(tmp) .* zero_e;

		tmp = e2 == 0; % cos(az) = 0
		zero_e = (-ry/2 <= us(tmp)) & (us(tmp) <= ry/2);
		ll(tmp) = ll(tmp) .* zero_e;

		tmp = e3 == 0; % sin(polar) = 0
		zero_e = (-rz/2 <= vs(tmp)) & (vs(tmp) <= rz/2);
		ll(tmp) = ll(tmp) .* zero_e;

if 0 % old way.  todo cut.  and verify e1 and e2 cases above
		if (e3 == 0) % sin(polar) = 0; so line along z
			zero_e = (-rz/2 <= vs) & (vs <= rz/2);
			ll = ll .* zero_e;
		end

		if (e1 == 0) % sin(az) = 0;
			zero_e = (-rx/2 <= us) & (us <= rx/2);
			ll = ll .* zero_e;
		end

		if (e2 == 0) % cos(az) = 0;
			zero_e = (-ry/2 <= us) & (us <= ry/2);
			ll = ll .* zero_e;
		end
end

		proj(:,:,ib) = proj(:,:,ib) + val * ll;
		clear proj_i
	end
end

end % cuboid_proj_do()


% cuboid_proj_test()
% internal test routine
function cuboid_proj_test

if 1
	params = [30 20 5  200 100 20  45 0 10];
	fun_proj = @(cg) cuboid_proj(cg, params, 'oversample', 2); % analytical
	fun_im = @(ig) cuboid_im(ig, params, 'oversample', 2); %, 'checkfov', true);

	ir_proj3_compare1(fun_proj, fun_im, 'chat', 1);
return % todo
end

down = 30;
cg = ct_geom('fan', 'ns', round(888/down), 'nt', 64, ...
	'na', 18, ...
	'ds', 1.0*down, 'dt', 1.1, ...
	'down', 1, ... % only downsample s and beta
	'dsd', 949, 'dod', 408, 'dfs', 0, ... % 3rd gen CT
	'dsd', inf, ... % parallel-beam
	'pitch', 0.5, ...
	'offset_t', 0.0, ...
	'offset_s', 0.25); % quarter detector
%	'dsd', 949, 'dod', 408, 'dfs', inf, ... % flat detector

ig = image_geom('nx', 512/2^2, 'nz', 64, 'fov', 512, 'dz', 0.625);
%cg.plot3(ig);

params = [30 20 5  200 100 20  45 0 10];
x = cuboid_im(ig, params, 'oversample', 2);

im plc 2 2
im(1, x), cbar

ya = cuboid_proj(cg, params, 'oversample', 2);
im(2, cg.s, cg.t, ya), cbar
titlef('analytical cone-beam projections, dfs=%g', cg.dfs)
xlabel s, ylabel t

%im clf, im(cg.s, cg.ad, permute(ya, [1 3 2])), cbar
%xlabel s, ylabel '\beta'

A = Gcone(cg, ig);
yd = A * x;

im(3, cg.s, cg.t, yd, 'discrete projections'), cbar

im(4, cg.s, cg.t, ya-yd), cbar
titlef('analytical - discrete')

max_percent_diff(ya, yd)

end % cuboid_proj_test()
