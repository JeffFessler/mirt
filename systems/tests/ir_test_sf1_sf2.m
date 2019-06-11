% ir_test_sf1_sf2.m
%
% Compare versions of separable footprint (sf1, sf2, sf3)
% for the case of large voxels and small rays.
%
% In this case, SF1 (ray-dependent scaling) is better than SF2,
% except in one region!?
%
% Copyright 2017-02-01, Jeff Fessler, University of Michigan


%% small system
if ~isvar('cg'), printm 'setup small'
	f.down_ig = 16;
	f.down_cg = 1;

	cg = ct_geom('ge2', ...
		'nt', 3*f.down_cg, ... % one row would suffice
	...%	'source_z0', -20, 'pitch', 0.3, ... % test helix
		'source_z0', 0, ...
		'pitch', 0, ...
		'orbit_start', 42, ...
		'na', 1*f.down_cg, ... % one view
...%		'dsd', inf, ... % par
...%		'dfs', inf, ... % flat
...%		'dsd', 3000, ... % nearly parallel
		'down', f.down_cg);

	ig = image_geom('nx', 512, 'ny', 480, 'nz', 528, ...
		'dx', 500/512, 'dy', '-dx', 'dz', 0.5, ...
		'dy', -0.4, ... % stress test!
		'nz', 1 * f.down_ig, ...
		'offset_x', 2.9, 'offset_y', 3.5, 'offset_z', -3.4, ... % stress
		'down', f.down_ig);

	tmp = logical(ig.zeros);
%	tmp(13:27,13:22,4:19) = true;
%	tmp(8:27,16,4:19) = true; % one 'iy' strip for testing
%	tmp(8:29,18,:) = true; % one 'iy' strip for testing
%	tmp(18,8:29,:) = true; % one 'ix' strip for testing
%	ig.mask = ig.unitv('c', [8 7 3]) > 0;
	tmp(6:end-4,6:end-4) = true; % big square
	ig.mask = tmp;

	if im, cg.plot(ig); end
end


%% analytical projection
if ~isvar('ya'), printm 'ya'
	x0 = single(ig.mask);
	cfun = @(c) mean(c(ig.mask)); % rect center
	wfun = @(c,d) abs(diff(minmax(c(ig.mask))) + abs(d)); % rect width
	cub = [	cfun(ig.xg) ...
		cfun(ig.yg) ...
		cfun(ig.zg) ...
		wfun(ig.xg, ig.dx) ...
		wfun(ig.yg, ig.dy) ...
		wfun(ig.zg, ig.dz) ...
		0 0 1];
	xa = cuboid_im(ig, cub);
	xa = round(xa); % trick: 0,1
	jf_equal(xa, x0)

%	im(x0)
%	im_toggle(x0, xa)
%	im(x0-xa)

	ya = cuboid_proj(cg, cub, 'oversample', 8);
end


%% examine trapezoids - todo: needs more work
if 0
	im(ig.x, ig.y, ig.mask)
	xc = ig.xg(ig.mask);
	yc = ig.yg(ig.mask);
	ar = cg.ar;
	hold on
	plot(xc, yc, 'o')
	hold off
	cos_a = cos(ar);
	sin_a = sin(ar);
	tproj0 = xc * cos_a + yc * sin_a;
	dsxy0 = cg.dso + xc * sin_a - yc * cos_a;
	tproj0_dsxy0 = sqrt(tproj0 .* tproj0 + dsxy0 .* dsxy0);
	dx_sin_2 = ig.dx * sin_a / 2.;
	dx_cos_2 = ig.dx * cos_a / 2.;
	dy_sin_2 = ig.dy * sin_a / 2.;
	n_dy_cos_2 = -ig.dy * cos_a / 2.;
	dsd_ds = cg.dsd / cg.ds;

	ratio = @(a,b) (tproj0 + a * dx_cos_2 + b * dy_sin_2) ...
			./ (dsxy0 + a * dx_sin_2 + b * n_dy_cos_2);
	if isinf(cg.dfs) % flat
		taufun = @(a,b) dsd_ds * ratio(a,b)';
	elseif ~isinf(cg.dsd) % arc
		taufun = @(a,b) dsd_ds * atan(ratio(a,b))';
	else
		fail 'not done'
	end

	taus = [taufun(-1,-1); taufun(+1,-1); taufun(-1,+1); taufun(+1,+1)];
	taus = sort(taus, 1);
	plot(taus')
	nr = 2001;
	rr = linspace(-300, 200, nr)';
	traps = zeros(nr, numel(xc));
	trapfun = @(taus, r) ...
		0 .* (r < taus(1)) + ...
		(r - taus(1)) ./ (taus(2) - taus(1)) ...
			.* (taus(1) < r & r <= taus(2)) + ...
		1 .* (taus(2) < r & r <= taus(3)) + ...
		(taus(4) - r) ./ (taus(4) - taus(3)) ...
			.* (taus(3) < r & r <= taus(4));
		
	for ix = 1:numel(xc)
		traps(:,ix) = trapfun(taus(:,ix), rr);
	end
	trapsum = sum(traps, 2);

	SFAmpRect2 = @(dx, dy, cos_phi, sin_phi) ...
		abs(dx * dy) ./ max(abs(dx * cos_phi), abs(dy * sin_phi));
	SFAmpRect = @(dx, dy, phi) SFAmpRect2(dx, dy, cos(phi), sin(phi));

	s_r = rr*cg.ds;
	phi_r = cg.gamma_s(s_r) + ar;
	dr = rad2deg(phi_r);
	da = rad2deg(cg.gamma + ar);
	scale_phi = SFAmpRect(ig.dx, ig.dy, phi_r);
	scale = scale_phi ./ sqrt(1 + s_r.^2 / cg.dsd.^2);
	plot(da, ya, '-', ...
		dr, ig.dx * trapsum, '-', ...
		dr, scale .* trapsum, '-o', ...
		dr, ig.dx*traps, '-')
	xlim([30,60])
return
end


if ~isvar('A3'), printm 'A'
	args = {cg, ig};
	A0 = Gcone(args{:}, 'type', 'rf1');
	A1 = Gcone(args{:}, 'type', 'sf1');
	A2 = Gcone(args{:}, 'type', 'sf2');
	A3 = Gcone(args{:}, 'type', 'sf3');
end



if 1
	cpu etic
	y0 = A0 * x0;
	cpu etoc rf1
	cpu etic
	y1 = A1 * x0;
	cpu etoc sf1
	y2 = A2 * x0;
	y3 = A3 * x0;

%	is = 240:260;
%	is = 150:320;
	is = 1:cg.ns;
	ifun = @(y) im(y(is,:));
	efun = @(y,s) im(y(is,:) - ya(is,:), [-1 1]*0.2, s);
	im plc 2 5
	ifun(ya), axis normal
	ifun(y0), axis normal
	ifun(y1), axis normal
	ifun(y2), axis normal
	ifun(y3), axis normal
	im subplot 7
	efun(y1, 'y0 - ya'), cbar h, axis normal
	efun(y1, 'y1 - ya'), cbar h, axis normal
	efun(y2, 'y2 - ya'), cbar h, axis normal
	efun(y3, 'y3 - ya'), cbar h, axis normal
end


if 1
	prompt
%	it = cg.nt/2 + 3;
	it = floor((cg.nt+1)/2);
	ia = 1;
	subplot(313)
	clf
	if isinf(cg.dsd)
		x = cg.s(is);
	else
		x = rad2deg(cg.gamma + cg.ar(ia));
	end
	plot(	x, ya(is,it,ia), 'k.-', ...
		x, y0(is,it,ia), 'md', ...
		x, y1(is,it,ia), 'r+', ...
		x, y2(is,it,ia), 'go', ...
		x, y3(is,it,ia), 'bx')
	xlabelf('$\phi$ [degrees]')
	ir_legend({'A', 'RF1', 'SF1', 'SF2', 'SF3'})
end
