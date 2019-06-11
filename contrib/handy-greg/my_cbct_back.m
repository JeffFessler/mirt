  function back = my_cbct_back(proj, cg, ig, varargin)
% function back = cbct_back(proj, cg, ig, varargin)
%|
%| cone-beam backprojector for feldkamp.m
%| in
%|	proj	[ns nt na]	cone-beam projection views
%| option
%|	'use_mex' 1|2|3		default: 1 mex with loop in mex
%|					2 mex with loop in matlab
%|					3 mex with "st" data order (slower)
%|					0 no mex; use matlab
%|	'ia_skip'		default: 1
%|	'offset_source'		default: 0
%|	'nthread'		default: jf('ncore')
%| out
%|	back	[nx ny nz]	back projection result
%|
%| Copyright 2004-8-28 Jeff Fessler, University of Michigan

if nargin == 1 && streq(proj, 'test'), cbct_back_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.ia_skip = 1;
arg.use_mex = 1;
arg.offset_source = 0;
arg.nthread = jf('ncore');
arg.back_call = @fdk_mex_call;
arg.scale_dang = true;
arg = vararg_pair(arg, varargin);

if arg.use_mex
	if cg.pitch ~= 0 || any(cg.source_zs ~= 0)
		fail('helix not yet supported')
	end
	back = cbct_back_mex(proj, cg.ns, cg.nt, cg.na, ...
		cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
		cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
		ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
		[ig.offset_x ig.offset_y ig.offset_z], ...
		arg.ia_skip, arg.scale_dang, ...
		arg.use_mex, arg.nthread, arg.back_call);

else

	back = cbct_back_mat(proj, cg, cg.ns, cg.nt, cg.na, ...
		cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
		cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
		cg.source_zs, ...
		ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
		[ig.offset_x ig.offset_y ig.offset_z], ...
		arg.ia_skip, arg.scale_dang);
end

end % cbct_back()


%
% cbct_back_mex()
% mex backprojector
%
function img = cbct_back_mex(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, scale_dang, ...
	use_mex, nthread, back_call)

[nx ny] = size(mask);
ia_list = 1:ia_skip:na;
betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles
proj = single(proj);
nthread = int32(nthread);

switch use_mex
case 1 % loop in mex
	proj = permute(proj, [2 1 3]); % sta -> tsa
	img = back_call('fdk,ts,back', ...
			int32([nx ny nz]), [dx dy dz], ...
			offset_xyz, uint8(mask), ...
			dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
			proj(:,:,ia_list), betas(ia_list), nthread);
%	img = double6(img);
	img = permute(img, [2 3 1]); % zxy -> xyz

case 2 % loop in matlab (for testing)
	proj = permute(proj, [2 1 3]); % sta -> tsa
	img = 0;
	for ia=ia_list
		ticker(mfilename, ia, na)
		% note: 2006-5-30: replaced -dy with dy
		tmp = back_call('fdk,ts,back', ...
			int32([nx ny nz]), [dx dy dz], ...
			offset_xyz, uint8(mask), ...
			dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
			proj(:,:,ia), betas(ia), nthread);
%		tmp = double6(tmp);
		img = img + tmp;
	end
	img = permute(img, [2 3 1]); % zxy -> xyz

case 3 % fdk,st (for testing only - slower!)
	warn 'fdk,st does not handle projection view edges completely' 
	img = back_call('fdk,st,back', ...
			int32([nx ny nz]), [dx dy dz], ...
			offset_xyz, uint8(mask), ...
			dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
			proj(:,:,ia_list), betas(ia_list), nthread);
%	img = double6(img);

otherwise
	fail 'bug'
end

if scale_dang % final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
end

end % cbct_back_mex()


%
% cbct_back_mat()
% matlab backprojector (slower)
%
function img = cbct_back_mat(proj, cg, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	source_zs, ...
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, scale_dang)

if any(source_zs ~= 0)
	warn('helix not yet tested')
end

[nx ny] = size(mask);
betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles

% precompute as much as possible
wx = (nx-1)/2 + offset_xyz(1);
wy = (ny-1)/2 + offset_xyz(2);
wz = (nz-1)/2 + offset_xyz(3);
[xc yc] = ndgrid(([0:nx-1] - wx) * dx, ([0:ny-1] - wy) * dy);
zc = ([0:nz-1] - wz) * dz;

if 0 % limit back-projection to FOV?  removed 2008-10-9
	rr = sqrt(xc.^2 + yc.^2); % [nx,ny]
	smax = ((ns-1)/2-abs(offset_s)) * ds; % maximum detector s coordinate

	if isinf(dfs)
		gamma_max = atan(smax/dsd);
	elseif dfs == 0
		gamma_max = smax / dsd;
	end

	rmax = dso * sin(gamma_max);
	mask = mask & (rr < rmax);
end
clear wx wy wz rr smax rmax

xc = xc(mask); % [np] pixels within mask
yc = yc(mask);

ws = (ns+1)/2 + offset_s; % trick: +1 because matlab starts from 1
wt = (nt+1)/2 + offset_t;

% loop over slices
img = zeros([size(mask) nz]);
sdim = [ns+3 nt+3]; % trick: extra zeros saves indexing in loop
proj1 = zeros(sdim);
ticker reset

% Greg Handy stuff:
num_turns = orbit/360;
halfNumAngles = round((na/num_turns)/2);
myPitch = cg.pitch * cg.nt * dso / dsd * cg.dt;
h = myPitch / (2*pi);
dist = .5 * myPitch;
zindex = 1;

for iz=1:nz

	% Greg Handy stuff:
	% Just trying to limit the number of angles used
	lower_limit = zc(iz) - dist;

	% enter into the acceptable range for the z-slices
	while zindex < size(cg.source_zs,1)+1 && cg.source_zs(zindex) < lower_limit 
		zindex = zindex + 1;
	end
	lowerA = zindex - halfNumAngles;
	upperA = zindex + halfNumAngles;
	
	% a check to prevent error; images without enough views will be poor 
	if lowerA <= 0
		lowerA = 1;
	end
	if upperA > na
		upperA = na;
	end

	% loop over each projection angle
	img2 = 0;
	for ia=lowerA:ia_skip:upperA
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc * cos(beta) + yc * sin(beta);
		y_betas = dso - (-xc * sin(beta) + yc * cos(beta));

		% detector indices
		if isinf(dsd) || isinf(dso) % par
			mag = ones(size(y_betas));
		else
			mag = dsd ./ y_betas;
		end

		if isinf(dfs) ... % flat
			|| isinf(dsd) || isinf(dso) % par
			sprime = mag .* x_beta;
		elseif (dfs == 0) % arc
			r_loop = x_beta - offset_source;
			sprime = dsd * atan2(r_loop, y_betas);
		end

		tprime = mag * (zc(iz) - source_zs(ia));
		bs = sprime / ds + ws;
		bt = tprime / dt + wt;

		% bi-linear interpolation:
		is = floor(bs); % left bin
		it = floor(bt);

		wr = bs - is;	% left weight
		wl = 1 - wr;	% right weight
		wu = bt - it;	% upper weight
		wd = 1 - wu;	% lower weight

		ibad = (is < 0) | (is > ns) | (it < 0) | (it > nt);
		is(ibad) = ns+1; % trick! point at harmless zeros
		it(ibad) = nt+1;

		proj1(1+[1:ns],1+[1:nt]) = proj(:,:,ia); % trick: left side
		p1 =	wl .* proj1(sub2ind(sdim, is+1,it+1)) + ...
			wr .* proj1(sub2ind(sdim, is+2,it+1));
		p2 =	wl .* proj1(sub2ind(sdim, is+1,it+2)) + ...
			wr .* proj1(sub2ind(sdim, is+2,it+2));

		p0 = wu .* p1 + wd .* p2; % vertical interpolation

		if isinf(dfs) ... % flat
			|| isinf(dsd) || isinf(dso) % par
			p0 = p0 .* mag.^2; % back-projection weighting for flat
		elseif dfs == 0 % arc
			p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + y_betas.^2);
		end

		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask);
end % iz

if scale_dang % final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na*ia_skip)) * img;
end

end % cbct_back_mat()


%
% fdk_mex_call()
%
function out = fdk_mex_call(varargin)
if exist('fdk_mex') == 3
	printm 'using fdk_mex'
	out = fdk_mex(varargin{:});
elseif exist('jf_mex') == 3
	out = jf_mex(varargin{:});
else
	fail('bug: neither fdk_mex nor jf_mex found')
end
end % fdk_mex_call()


%
% cbct_back_test()
% compare various versions to check consistency
%
function cbct_back_test

down = 8; % fast test
% todo: test parallel, flat, arc
%cg = ct_geom('ge1', 'nt', 800, 'na', 50*down, 'dsd', inf, 'down', down);
cg = ct_geom('ge1', 'nt', 800, 'na', 50*down, 'down', down);
%cg.dt = -cg.dt; % todo: handle this case
ig = image_geom('nx', 512, 'ny', 496, 'nz', 480, 'fov', 500, ...
	'mask', 'all-but-edge', 'down', down);

ell = [ ... % somewhat realistic phantom object
	[30 10 10	150 150 280	0 0 1000]; % 30cm diam
	[80 10 10	50 50 30	0 0 300]; % bone-like inserts
	[-10 -40 75	40 40 40	0 0 600];
	[-10 80 -20	30 30 30	0 0 600];
];

if 1
	proj = ellipsoid_proj(cg, ell);
	proj = fdk_filter(proj, 'ramp', cg.dsd, cg.dfs, cg.ds);
else
	proj = cg.zeros;
	proj(:,:,9) = 1;
end
% im clf, im(proj, 'true projections'), cbar, return

if 0 % zero outer edges of projection for comparing ,st to ,ts
	proj([1 cg.ns], :, :) = 0;
	proj(:, [1 cg.nt], :) = 0;
	f.compare_st = true;
else
	f.compare_st = false;
end

	ia_skip = 3; % stress test and makes it faster too
	args = {proj, cg, ig, 'ia_skip', ia_skip, 'scale_dang', 0};
if 1
	back0 = cbct_back(args{:}, 'use_mex', 0);
%	im(back0), return
end

if exist('jf_mex') == 3
	for ii = 1:3
		use_mex = ii;
		printm('testing jf_mex with use_mex=%d', use_mex)
		back1{ii} = cbct_back(args{:}, ...
			'use_mex', use_mex, 'back_call', @jf_mex);
	end
	printm('mpd mex1 vs mat: %g%%', max_percent_diff(back1{1}, back0))
	printm('mpd mex1 vs mex2: %g%%', max_percent_diff(back1{1}, back1{2}))
	if f.compare_st % only if edges are zeroed
	printm('mpd mex1 vs mex3: %g%%', max_percent_diff(back1{1}, back1{3}))
	end
end

if exist('fdk_mex') == 3
	printm('found fdk_mex:\n %s', which('fdk_mex'))
	back2 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @fdk_mex);
	printm('mpd jf_mex vs fdk_mex: %g%%', max_percent_diff(back1{1}, back2))

	if 1
		im plc 1 3
		iz = 1:ig.nz;
		im(1, back1{1}(:,:,iz))
		im(2, back2(:,:,iz))
		err = back2(:,:,iz) - back1{1}(:,:,iz);
		im(3, err)
	end
end

end % cbct_back_test()
