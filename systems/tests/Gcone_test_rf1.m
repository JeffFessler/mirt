% Gcone_test_rf1.m
% test the Gcone object for rf1
% Copyright 2017-01-29, Jeff Fessler, University of Michigan

% todo: first compare for square pixels vs DD
% todo: then compare for rect pixels vs sf2

if exist('dd_ge1_mex') ~= 3 % UM only
	fail 'need DD'
end

f.class = 'fatrix2';
f.is_zxy = false;
%f.is_zxy = true; % todo
if f.is_zxy
	to_zxy = @(x) permute(x, [3 1 2]);
else
	to_zxy = @(x) x;
end


%% small system
if 1 || ~isvar('Asf2'), printm 'setup small'
	f.down_ig = 16;
	f.down_cg = 1;

	cg = ct_geom('ge2', 'nt', 64, ...
		'nt', 3, ...
		'source_z0', -20, 'pitch', 0.5, ... % test helix
		'orbit_start', 40, 'na', 1, ... % test 1 view
		'down', f.down_cg);
%'dfs', inf, ... % flat
%'dsd', inf, ... % parallel

	ig = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 500/512, 'dy', '-dx', 'dz', 0.5, ...
		'offset_x', 2.9, 'offset_y', 3.5, ...
		'offset_z', -3.4, ...
		'mask', 'all-but-edge', ... % trick: for hct2
		'down', f.down_ig);
%	ig.mask = ig.circ > 0;
	tmp = logical(ig.zeros);
%	tmp(13:27,13:22,4:19) = true;
	tmp(8:27,16,4:19) = true; % one 'iy' strip for testing
	ig.mask = tmp;
%	im(ig.mask), return

	args = {cg, ig, 'type', 'sf2', 'class', f.class};
	args = {args{:}, 'zxy', f.is_zxy};
	Asf2 = Gcone(args{:});
end


if 1, printm 'rf1' % test rf1 vs dd2
	Arf1 = Gcone(args{:}, 'type', 'rf1');
	if isinf(cg.dsd) % parallel
		Add2 = Gcone(args{:}, 'type', 'sf1');
	else
		Add2 = Gcone(args{:}, 'type', 'dd2');
	end

	x0 = single(ig.mask);
%	x0 = ig.unitv('c', [ 8 6 4 ]);
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
%	im(x0)
%	im_toggle(x0, xa)
%	im(x0-xa)

	ya = cuboid_proj(cg, cub, 'oversample', 4);

	if 1 % back
		xrf1 = Arf1' * ya;
		xdd2 = Add2' * ya;
		im plc 1 3
		im(xrf1), cbar
		im(xdd2), cbar
		im(xrf1 - xdd2), cbar
		max_percent_diff(xrf1, xdd2)
		prompt
	end

	if 1
		cpu etic
		yd = Add2 * x0;
		cpu etoc dd
	else
	%	yd = cg.zeros;
		yd = ya;
	end

	cpu etic
	y2 = Asf2 * x0;
	cpu etoc sf2

	if 1
		cpu etic
		y0 = Arf1 * x0; % todo: crashes matlab
		cpu etoc rf1
	end

	is = 1:cg.ns;
	efun = @(y,s) im(y(is,:) - ya(is,:), [-1 1]*0.2, s);
	im plc 2 4
	im(ya), axis normal
	im(yd), axis normal
	im(y0), axis normal
	im(y2), axis normal
	im subplot 5
	cg.plot(ig);
	im subplot 6
	efun(yd, 'yd-ya'), cbar h, axis normal % todo: unexpected large difference
	efun(y0, 'y0-ya'), cbar h, axis normal % small (1%) as expected
	efun(y2, 'y2-ya'), cbar h, axis normal

	max_percent_diff(y0, yd)
end


if 1
	prompt
	ii = cg.s;
	it = floor((cg.nt+1)/2);
	ia = 1;
	clf, plot( ...
		ii, ya(:,it,ia), 'b.-', ...
		ii, yd(:,it,ia), 'r.-', ...
		ii, y0(:,it,ia), 'md', ...
		ii, y2(:,it,ia), 'k.-')
	ir_legend({ 'Analytical', Add2.type, 'RF1', 'SF2'})
end
