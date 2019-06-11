% Gtomo2_table_test
% Test the Gtomo2_table obect by comparing with exact methods
% Jeff Fessler

% paces and adjoint
if ~isvar('At2'), printm 'test Gtomo2_table and its adjoint'
	sgt = sino_geom('par', 'nb', 20, 'na', 10);
	igt = image_geom('nx', 7, 'ny', 12, 'dx', 1);
	igt.mask = igt.circ > 0;
	table = {'square/strip', 'chat', 0, 'Ltab', 1000, ...
		'strip_width', sgt.dr};
	At1 = Gtomo2_table(sgt, igt, table, 'nthread', 1);
	At2 = Gtomo2_table(sgt, igt, table);
	tester_tomo2(At1, igt.mask, 'A2', At2)
	test_adjoint(At1);
	test_adjoint(At2);
%	clear sgt igt At1 At2 table
end

redo = 0;
% analytical sino
if redo || ~isvar('ya'), printm 'ya'
	down = 8;
%	down = 2; % for threading check
	ig = image_geom('nx', 512, 'ny', 496, 'fov', 500, 'down', down);
	ig.mask = ig.circ > 0;

	sg = sino_geom('par', 'nb', 888, 'na', 984, ...
		'dr', 1.0, 'offset_r', 0.25, ...
		'strip_width', 'd', ...
		'orbit', 0*360+90, 'orbit_start', -15*0, 'down', down);

	f.mask = [test_dir 'mask.fld'];
	fld_write(f.mask, ig.mask, 'type', 'byte');

	[xa ell] = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	im clf, im(xa, 'xa')

	ya = ellipse_sino(sg, ell, 'oversample', 8);
	im(ya, 'ya'), cbar
prompt
end


% Atab
if redo || ~isvar('Atab'), printm 'setup Atab'
	if 0 % linear interpolation (cf system 9)
		if sg.strip_width ~= 0, fail 'change sg', end
		f.system = 9;
		table = {'linear', 'chat', 0, 'Ltab', 1000};

	elseif 1 % square-pixel / strip-integral
		f.system = 2;

		f.table = {'square/strip', 'chat', 0, 'Ltab', 1000, ...
			'strip_width', sg.dr};
	end

	cpu etic
	Atab = Gtomo2_table(sg, ig, f.table, 'nthread', 1);
	cpu etoc 'Atab precompute time:'
	Atab2 = Gtomo2_table(sg, ig, f.table);
prompt
end

if 0, printm 'threading check'
	cpu etic
	y1 = Atab * xa;
	cpu etoc 'Atab1 *:'

	cpu etic
	y2 = Atab2 * xa;
	cpu etoc 'Atab2 *:'
	jf_equal(y1, y2)
end

if exist('dd_ge1_mex') == 3 % DD approx, UM only
	if ~isvar('As')
		t = f.table; t{5} = 200;
		As = Gtomo2_table(sg, ig, t); % exact
		Ad = Gtomo2_table(sg, ig, {'dd2', t{2:end}});
	end
	A45 = Gtomo2_table(sg, ig, {'la45', t{2:end}});
	K = As.arg.tab_opt.Ktab;
	L = As.arg.tab_opt.Ltab;
	[kk ll] = ndgrid(0:K-1, 0:L-1);
	fs = double6(As.arg.table);
	fd = double6(Ad.arg.table);
	f4 = double6(A45.arg.table);
	minmax(sum(fs, 1))
	minmax(sum(fd, 1))
	minmax(sum(f4, 1))
	fs = reshape(fs, L*K, []);
	fd = reshape(fd, L*K, []);
	f4 = reshape(f4, L*K, []);
	ti = col(As.arg.t_fun(kk, ll));
	[ti ii] = sort(ti);
	fs = fs(ii, :);
	fd = fd(ii, :);
	f4 = f4(ii, :);
	fde = fd - fs;
	f4e = f4 - fs;
	if im
		im plc 2 2
		alist = [21 32 45];
		for ii=1:3
			ia = imin(abs(sg.ad - alist(ii)));
			scale = 1;% / max(fs(:));
			im('subplot', ii)
			plot(	ti, scale * fd(:,ia)-0*fs(:,ia), 'c--', ...
				ti, scale * fs(:,ia)-0*fs(:,ia), 'r-', ...
				ti, scale * f4(:,ia)-0*fs(:,ia), 'y:')
			legend('DD', 'Exact', 'LA45')
			xlabel 'radial position'
			ylabel 'footprint'
			titlef('projection view angle %g', sg.ad(ia))
		end

%		im([1*fs, 1*fd, fde]), cbar
%		im pl 2 1,
		im subplot 4
		plot(	sg.ad, max(abs(fde), [], 1), 'c-', ...
			sg.ad, mean(abs(fde), 1), 'c-o', ...
			sg.ad, mean(abs(f4e), 1), 'y-o', ...
			sg.ad, max(abs(f4e), [], 1), 'y-')
		axis([0 180 0 1])
		legend('DD', 'DD err', 'new err', 'new')
		xtick([0:45:360])
		xlabel 'projection angle'
		ylabel 'projection error'
		titlef('error for dx=%g, dr=%g', ig.dx, sg.dr)
	end
prompt
end

% dsc
if redo || ~isvar('Adsc'), printm 'setup Adsc'
	Adsc = Gtomo2_dscmex(sg, ig);
	if f.system == 9 % untested
		scale = sg.orbit / 180 / (2*pi) * (sg.dr)^2 / ig.dx; % trick:
		Adsc = scale * Adsc;
	end
prompt
end

% wtf
if redo || ~isvar('Awtc'), printm 'setup Awtc/Awtr'
	if ig.nx*ig.ny <= 2^16 && sg.nb*sg.na <= 2^16
		Awtr = Gtomo2_wtmex(sg, ig);
	end
	if ig.nx*ig.ny <= 2^16 && sg.nb*sg.na <= 2^16
		Awtc = Gtomo2_wtmex(sg, ig, 'grouped', 'col');
	end
end


if 1, printm 'forward'
	if 0
%		x = double(ig.mask);
		x = zeros(size(ig.mask));
		x(round(ig.nx/4), round(ig.ny/3)) = 1;
	else
		x = xa;
	end
	yt = Atab * x;
	yd = Adsc * x;
	yc = Awtc * x;
	yr = Awtr * x;

	cpu etic
	yc = Awtc * x;
	cpu etoc 'Awtc *:'

	cpu etic
	yr = Awtr * x;
	cpu etoc 'Awtr *:'

	cpu etic
	yd = Adsc * x;
	cpu etoc 'Adsc *:'

%profile on
	cpu etic
	yt = Atab * x;
	cpu etoc 'Atab1 *:'
%profile report

	cpu etic
	yt2 = Atab2 * x;
	cpu etoc 'Atab2 *:'

	if max_percent_diff(yt, yt2), error 'thread bug', end
	max_percent_diff(yt, yd)
	max_percent_diff(yt, yr)
	max_percent_diff(yr, yc)
	im clf, im([yt; yd])
%	plot(1:na, sum(yt,1), '-', 1:na, sum(yd,1), '.')
	im(yd-yt), cbar
%	plot(1:nb, yd(:,1), '.', 1:nb, yt(:,1), 'o')
%	sum(yt(:)) ./ sum(yd(:))
%	worst = imax(abs(yt-yd),2)
%	rad2deg(ang1(worst(2)))
prompt
end


% test back projection
if 1, printm 'back'
	y = ya;
%	y = ones(size(ya));
%	y = zeros(size(ya));
%	y(round(nb/2), :) = 1;
	xd = Adsc' * y;
	xc = Awtc' * y;
	xr = Awtr' * y;
	xt = Atab' * y;
	im([xt; xr; xd])

	cpu etic
	xc = Awtc' * y(:);
	cpu etoc 'Awtc back:'
	xc = ig.embed(xc);

	cpu etic
	xr = Awtr' * y(:);
	cpu etoc 'Awtr back:'
	xr = ig.embed(xr);

	cpu etic
	xd = Adsc' * y;
	cpu etoc 'Adsc back:'

	cpu etic
	xt = Atab' * y;
	cpu etoc 'Atab1 back:'

	cpu etic
	xt2 = Atab2' * y;
	cpu etoc 'Atab2 back:'
	max_percent_diff(xt, xt2)

	max_percent_diff(xt, xd)
	max_percent_diff(xt, xr)
	max_percent_diff(xc, xr)
prompt
end

% hereafter requires DD
if 2 ~= exist('Gtomo_dd') || exist('dd_ge1_mex') ~= 3
	return
end

if redo || ~isvar('Add'), printm 'setup Add'
	cpu etic
	Add = Gtomo_dd(sg, ig);
	cpu etoc 'Ad precompute time:'
end

if 1, printm 'tab vs dd'
	% cache warm up
	yt = Atab * xa;
	yd = Add * xa;
	cpu etic
	yt = Atab * xa;
	cpu etoc 'tab'
	cpu etic
	yd = Add * xa;
	cpu etoc 'dd'
	max_percent_diff(ya, yt)
	max_percent_diff(ya, yd)
	max_percent_diff(yt, yd)
return
end

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0 && ~isvar('An'), printm 'setup An'
warning 'todo: update nufft'
	An = Gtomo_nufft(ig.mask, [nb na], ...
		'chat', 1, ...
		fan_arg{:}, ...
		'xscale', f.xscale, 'yscale', f.yscale, ...
		'orbit', f.orbit, ...
		'orbit_start', f.orbit_start, ...
		'pixel_size', f.pixel_size, ...
		'ray_spacing', f.ray_spacing, ...
		'strip_width', f.ray_spacing, ...
		'interp', {'table', 2^11, 'minmax:kb'});
	printf('An precompute time = %g', toc)
end
