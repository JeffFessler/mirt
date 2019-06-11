  function proj = ellipsoid_proj(cg, params, varargin)
%|function proj = ellipsoid_proj(cg, params, varargin)
%|
%| Compute set of 2d line-integral projection views of ellipsoid(s).
%| Works for both parallel-beam and cone-beam geometry.
%|
%| in
%|	cg			ct_geom()
%|	params [ne 9]		ellipsoid parameters:
%|			[x_center y_center z_center  x_radius y_radius z_radius
%|				xy_angle_degrees z_angle_degrees  amplitude]
%| options
%|	oversample		over-sampling factor for emulating "strips"
%|				(to account for finite detector size)
%|
%| out
%|	proj	[ns nt na]	projection views
%|
%| Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test'), ellipsoid_proj_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = ellipsoid_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % ellipsoid_proj()


% ellipsoid_proj_do()
function proj = ellipsoid_proj_do(params, ss, tt, ...
		beta, ... % [radians]
		source_zs, dso, dod, dfs, oversample)

if size(params, 2) ~= 9, error '9 parameters per ellipsoid', end

if oversample > 1
	ds = ss(2) - ss(1);
	dt = tt(2) - tt(1);
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6)
		error 'uniform spacing required for oversampling'
	end
	No = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No-1):2:(No-1)]'/(2*No)*ds, ss(:)'); % [No ns]
	tt = outer_sum([-(No-1):2:(No-1)]'/(2*No)*dt, tt(:)'); % [No nt]
	proj = ellipsoid_proj_do(params, ss(:), tt(:), beta, source_zs, dso, dod, dfs, 1);
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
	azim0 = 0;
	polar = 0;

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

% loop over ellipsoids
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
%		az = beta(ib) + azim0 - xang; % assume source rotate in xy plane
		az = beta(ib) + azim0; % correction due to Lei Zhu of Stanford

		% shift property of 3D transform:
		cz_eff = cz - source_zs(ib); % center relative to source
		ushift = cx * cos(az) + cy * sin(az);
		vshift = (cx * sin(az) - cy * cos(az)) .* spolar + cz_eff * cpolar;

		az = az - xang; % correction due to Lei Zhu of Stanford
		p1 = (uu-ushift) .* cos(az) + (vv-vshift) .* sin(az) .* spolar;
		p2 = (uu-ushift) .* sin(az) - (vv-vshift) .* cos(az) .* spolar;
		p3 = (vv-vshift) .* cpolar;

		e1 = -sin(az) .* cpolar;
		e2 = cos(az) .* cpolar;
		e3 = spolar;

		A = e1.^2 / rx^2 + e2.^2 / ry^2 + e3.^2 / rz^2;
		B = p1 .* e1 / rx^2 + p2 .* e2 / ry^2 + p3 .* e3 / rz^2;
		C = p1.^2 / rx^2 + p2.^2 / ry^2 + p3.^2 / rz^2 - 1;

		proj(:,:,ib) = proj(:,:,ib) + 2 * val * sqrt(B.^2 - A.*C) ./ A;
	end
end

% trick: anywhere proj of a single ellipsoid is imaginary, the real part is 0.
proj = real(proj);
end % ellipsoid_proj_do()


% ellipsoid_proj_test()
function ellipsoid_proj_test

ell = [20 0*50 -40 200 100 50 90 0 10;
        0 50 100 80 80 20 0 0 10];

fun_proj = @(cg) ellipsoid_proj(cg, ell, 'oversample', 2); % analytical
fun_im = @(ig) ellipsoid_im(ig, ell, 'oversample', 2, 'checkfov', true);

ir_proj3_compare1(fun_proj, fun_im, 'chat', 1);

end % ellipsoid_proj_test()
