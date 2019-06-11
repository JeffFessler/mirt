% Gtomo2_moj_test2
% Test the Mojette case of the Gtomo2_table object
% for fan-beam forward projection
% Copyright 2005-12-14, Jeff Fessler, University of Michigan

redo = 0;
im pl 2 3
% analytical sino
if redo || ~isvar('ya'), printm 'ya'
	down = 2;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	sgf = sino_geom('fan', 'nb', 888, 'na', 984/8, ...
		'ds', 1.02, 'offset_s', 0.25, ...
		'dsd', 949, 'dod', 408, ...
		'orbit', 360, 'orbit_start', -15, 'down', down);

	ig.mask = ellipse_im(ig, [0 0 250*[1 1] 0 1], 'oversample', 2) > 0;
	f.mask = [test_dir 'mask.fld'];
	fld_write(f.mask, ig.mask, 'type', 'byte');

	ell = [20 10 200 150 10 1];
	[xa ell] = ellipse_im(ig, ell, 'oversample', 3);
	im plc 2 3, clim = [0.9 1.1];
	im(1, xa, 'xa', clim), cbar

	ya = ellipse_sino(sgf, ell, 'oversample', 8);
	im(4, ya, 'ya'), cbar
prompt
end

% Gmpf (mojette: parallel to fan)
if redo || ~isvar('Gmpf'), printm 'setup Gmpf'
	f.table = {'mojette,square/strip'};

%	f.nr = ceil(sqrt(2)/2 * max(ig.nx,ig.ny))*2;
	f.nr = 4 + 2 * max(ig.nx,ig.ny);
	sgm = sino_geom('moj', 'nb', f.nr, 'na', sgf.na, ...
		'dx', ig.dx, 'offset_r', 0, ... % must be 0 for 180->360
		'orbit', 180, 'orbit_start', 0);

	cpu etic
	Gmoj = Gtomo2_table(sgm, ig, f.table);
	cpu etoc 'Gmoj precompute time:'

	Gpf = rebin_sino([], sgm, sgf, 'ob', 1);
	Gmpf = Gpf * Gmoj;
prompt
end

% dsc
if redo || ~isvar('Gdsc'), printm 'setup Gdsc'
	args = arg_pair('system', 14, 'nx', ig.nx, 'ny', ig.ny, ...
		'nb', sgf.nb, 'na', sgf.na, ...
...%		'support', 'all', ...
		'support', ['file ' f.mask], ...
		'pixel_size', ig.dx, ...
		'orbit', sgf.orbit, 'orbit_start', sgf.orbit_start, ...
		'src_det_dis', sgf.dsd, 'obj2det_x', sgf.dod, ...
		'ray_spacing', sgf.ds, ...
		'channel_offset', sgf.offset_s, ... % fix: -?
		'flip_y', 1, 'scale', 0);

	% todo:
	% Gdsc = Gtomo2_dscmex(sgf, igf, ...?);

	% system object
	Gdsc = Gtomo2_dscmex(args);
end

if redo || ~isvar('Gdd'), printm 'setup Gdd'
	cpu etic
	Gdd = Gtomo_dd(sgf, ig);
	cpu etoc 'Gd precompute time'
end

if redo || ~isvar('Gn'), printm 'setup Gn'
	cpu etic
	Gn = Gtomo_nufft(ig.mask, [sgf.nb sgf.na], ...
		'chat', 0, ...
		'dis_src_det', sgf.dsd, 'dis_iso_det', sgf.dod, ...
		'yscale', 1, ...
		'orbit', sgf.orbit, 'orbit_start', sgf.orbit_start, ...
		'dx', ig.dx, ...
		'ds', sgf.ds, 'offset_s', sgf.offset_s, 'strip_width', sgf.ds, ...
		'interp', {'table', 2^11, 'minmax:kb'});
	cpu etoc 'Gn precompute time'
end

if redo || ~isvar('Gw'), printm 'setup Gw'
	Gw = Gtomo2_wtmex(sgf, ig, 'grouped', 'col');
end

% test forward projection
if 1, printm 'proj'
	if 0 % cache warm up
		Gmpf(:,1);
		Gdd(:,1);
		Gn(:,1);
	end

	xr = xa;
	xr = ig.unitv('c', [20 5]);

	cpu etic
	ym = Gmpf * xr;
	cpu etoc 'mpf proj time'

	cpu etic
	yd = Gdd * xr;
	cpu etoc 'dd proj time'

	cpu etic
	yn = Gn * xr;
	cpu etoc 'nufft proj time'

	cpu etic
	yw = Gw * xr;
	cpu etoc 'wtmex proj time'

%	yr = ya;
	yr = yw;

	im(1, yn, 'yn'), cbar, cbar, axis normal
	im(2, ym, 'ym'), cbar, cbar, axis normal
	im(3, yd, 'yd'), cbar, cbar, axis normal
	im(4, yn-yr, 'yn-yr'), cbar, axis normal
	im(5, ym-yr, 'ym-yr'), cbar, axis normal
	im(6, yd-yr, 'yd-yr'), cbar, axis normal

	max_percent_diff(yr, ym)
	max_percent_diff(yr, yd)
	max_percent_diff(yr, yn)
%	max_percent_diff(ym, yd)
	nrms(ym, yr)
	nrms(yd, yr)
	nrms(yn, yr)
prompt
end

% test back projection
if 1, printm 'back'
	if 0 % cache warm up
		Gmpf(:,1);
		Gdd(:,1);
		Gn(:,1);
	end

	yy = zeros(size(ya)); yy(end/4,floor(sgf.na/4)) = 1; % one ray
	yy = ones(size(ya)); % sensitivity map

	cpu etic
%profile on
	xm = Gmpf' * yy;
%profile report
	cpu etoc 'mpf back time'

	cpu etic
	xd = Gdd' * yy;
	cpu etoc 'dd back time'

	cpu etic
	xn = Gn' * yy;
	cpu etoc 'nufft back time'

	if 0
		xn = xn - mean(xn(ig.mask)) * ig.mask;
		xm = xm - mean(xm(ig.mask)) * ig.mask;
		xd = xd - mean(xd(ig.mask)) * ig.mask;
	end
	im(1, xn, 'xn'), cbar
	im(2, xm, 'xm'), cbar
	im(3, xd, 'xd'), cbar
	im(4, xn-xm, 'xn-xm'), cbar
	im(5, xm-xd, 'xm-xd'), cbar
	im(6, xd-xn, 'xd-xn'), cbar

	max_percent_diff(xd, xm)
	max_percent_diff(xd, xn)
	nrms(xm, xd)
	nrms(xn, xd)
end

if 0, printm 'thread'
	tmp = sgm; tmp.na = 800;
	G1 = Gtomo2_table(tmp, ig, f.table, 'nthread', 1);
	G2 = Gtomo2_table(tmp, ig, f.table);

	arg = G1.arg;
	x = single(xr);
	cpu etic
%	ym = G1 * xr;
	y = jf_mex(arg.proj_str, arg.ptab_arg{:}, arg.table, ...
	        arg.mask8, arg.angles, x);
	cpu etoc 'moj1 proj time'

	arg = G2.arg;
	cpu etic
%	ym = G2 * xr;
	y = jf_mex(arg.proj_str, arg.ptab_arg{:}, arg.table, ...
	        arg.mask8, arg.angles, x);
	cpu etoc 'mojN proj time'
end
