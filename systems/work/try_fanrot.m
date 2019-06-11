% try_fanrot
% try fan-beam forward projection
% using hybrid rotation / precomputed 0 degree system

% this is on hold until rotmex bug is resolved
im(rotmex(ones(200, 'single'), 40.97561, int32(1)))
return


down = 4;
n.x = 512 / down;
n.y = n.x - 0;
n.b = 888 / down;
n.a = 984 / down;

f.orbit = 360;
f.orbit_start = 0;
f.orbit_start = 40.9756; f.orbit = 0.00002; % bug testing
f.dx = 500 / n.x; % 50cm fov
f.ds = 1.0 * down;

%
% ordinary completely dsc version
%
if ~isvar('G')
	f.dsd = 949;
	f.dod = 408;
	f.arg = arg_pair('system', 14, ...
		'nx', n.x, 'ny', n.y, 'nb', n.b, 'na', n.a, ...
		'pixel_size', f.dx, 'ray_spacing', f.ds, ...
		'orbit', f.orbit, 'orbit_start', f.orbit_start, ...
		'obj2det_x', f.dod, 'src_det_dis', f.dsd, ...
		'scale', 0);
	G = Gtomo2_dscmex(f.arg);
end

%
% dsc version for just view 0
%
if ~isvar('G0d')
	f.arg0 = arg_pair('system', 14, ...
		'nx', n.x, 'ny', n.y, 'nb', n.b, 'na', 1, ...
		'pixel_size', f.dx, 'ray_spacing', f.ds, ...
		'obj2det_x', f.dod, 'src_det_dis', f.dsd, ...
		'scale', 0);
	G0d = Gtomo2_dscmex(f.arg0);
end

%
% wtf version for just view 0
%
if ~isvar('G0w'), printm 'G0w'
	f.argw = sprintf(...
		'nx %d ny %d nb %d na 1 pixel_size %g ray_spacing %g', ...
		n.x, n.y, n.b, f.dx, f.ds);
	f.dis = sprintf('obj2det_x %g src_det_dis %g scale 0', f.dod, f.dsd);
	f.dsc = 't.dsc';
	f.wtf = 't.wtf';
	t = sprintf('wt -chat 0 dsc 14 %s %s >! %s', f.argw, f.dis, f.dsc);
	os_run(t)
        t = sprintf('echo y | wt gen %s', f.dsc);
        os_run(t)
	G0w = wtf_read(f.wtf);
end

if 1, printm 'rot0'
	filter = 3;
	chat = 0;
	cpu etic
	rotate_key = rotmex('init', int32(n.x), int32(n.y), int32(n.a), ...
		f.orbit, f.orbit_start, int32(filter), int32(chat));
	cpu etoc 'precompute rotate time'
end

%x = ones(n.x, n.y);
%x = zeros(n.x, n.y);
%x(round(n.x*2/3),end/2+1) = 1;
%x = double(G.mask);
x = shepplogan(n.x, n.y, 1);

if 0, printm 'G0d vs G0w'
	cpu etic
	y0d = G0d * x(:);
	t = cpu('etoc');
	printf('dsc 0 time %g * %d = %g', t, n.a, n.a*t)

	tic
	y0w = G0w * x(:);
	t = cpu('etoc');
	printf('wtf 0 time %g * %d = %g', t, n.a, n.a*t)
	plot([y0d y0w])

	max_percent_diff(y0d,y0w)
return
end

if 1, printm 'Gdsc * x'
	cpu etic
	y1 = reshape(G * x(:), n.b, n.a);
	cpu etoc 'all dsc time'
end

if 1, printm ''
%	profile off
%	profile on
	cpu etic
	y2 = proj_rotmex_wtf(x, G0w, n.a, rotate_key);
	cpu etoc 'rot time'
	max_percent_diff(y2,y1)
%	profile report
end

ang = (0:n.a-1)/n.a*f.orbit + f.orbit_start;
rad = 1:n.b;
im(121, rad, ang, y2), cbar
im(122, rad, ang, y1), cbar

if 0
%	profile off
%	profile on
	cpu etic
	y3 = proj_rot_wtf(x, G0w, n.a, filter);
	cpu etoc 'rot time'
	max_percent_diff(y3,y1)
	max_percent_diff(y3,y2)
%	profile report
end

% rotmex('free', key)
