 function back = cbct_back_helix(proj, cg, ig, varargin)
%function back = cbct_back_helix(proj, cg, ig, varargin)
%|
%| cone-beam backprojector for feldkamp.m
%|
%| in
%|	proj	[ns nt na]	cone-beam projection views
%|	cg	strum		ct_geom
%|	ig	strum		image_geom
%|
%| option
%|	'use_mex' 0|1|2|3	1 mex with loop in mex (default)
%|				2 mex with loop in matlab
%|				3 mex with "st" data order (slower)
%|				0 no mex; use matlab
%|	'ia_skip'		default 1
%|	'offset_source'		default 0
%|	'nthread'		default: jf('ncore')
%|	'back_call'		default @fdk_mex_call
%|	'scale_dang'		scale by "d angle"? default 1
%|	'extrapolate_t' 0|?	default 0
%|				(recommend 1.3 * nt/2 for ge1 geometry)
%|	'cone_par'		default 0 by DK
%|	'w3d'			default 0 by DK
%|
%| out
%|	back	[nx ny nz]	back projection result
%|
%| Copyright 2004-8-28 Jeff Fessler, University of Michigan
%| 2013-04-07 fixed wd,wu and improved extrapolation by Rebecca Malinas
%| 2014-07-16 JF added row extrapolation by padding before calling mex version
%|
%| 2014-11-10 helical version of cbct_back.m with view weighting
%| proposed in @tang:05:atd and @tang:06:atd, contibuted by Donghwan Kim

if nargin == 1 && streq(proj, 'test'), cbct_back_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.use_mex = []; % see below
arg.ia_skip = 1;
arg.offset_source = 0;
arg.nthread = jf('ncore');
arg.back_call = @fdk_mex_call;
arg.scale_dang = true;
arg.extrapolate_t = 0;
arg.cone_par = 0; % DK
arg.w3d = 0; % DK
arg = vararg_pair(arg, varargin);

if isempty(arg.use_mex)
	arg.use_mex = has_mex_jf; % default: usual mex version "1" if available 
end

if arg.use_mex
	back = cbct_back_mex(proj, cg.ns, cg.nt, cg.na, ...
		cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
		cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
		cg.source_zs, ... % DK
		cg.source_dz_per_view, cg.pitch, ... % DK
		ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
		[ig.offset_x ig.offset_y ig.offset_z], ...
		arg.ia_skip, arg.scale_dang, arg.extrapolate_t, ...
		arg.cone_par, arg.w3d, ... % DK
		arg.use_mex, arg.nthread, arg.back_call);

else
	back = cbct_back_mat(proj, cg.ns, cg.nt, cg.na, ...
		cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
		cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
		cg.source_zs, ...
		cg.source_dz_per_view, cg.pitch, ... % DK
		ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
		[ig.offset_x ig.offset_y ig.offset_z], ...
		arg.ia_skip, arg.scale_dang, arg.extrapolate_t, ...
		arg.cone_par, arg.w3d); %DK
end

end % cbct_back()


% cbct_back_mex()
% mex backprojector
function img = cbct_back_mex(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	source_zs, ... % DK
	source_dz_per_view, pitch, ... % DK
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, scale_dang, ...
	extrapolate_t, ...
	cone_par, w3d, ... % DK
	use_mex, nthread, back_call)

[nx ny] = size(mask);
ia_list = 1:ia_skip:na;
betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles
proj = single(proj);
nthread = int32(nthread);
na1_half = int32(ceil((na/ia_skip) / abs(orbit) * 180)); % DK
source_dz_per_rad = double(source_dz_per_view / (deg2rad(abs(orbit)) / (na/ia_skip))); % DK: H/(2pi)

if extrapolate_t > 0
	tmp1 = repmat(proj(:, 1, :), [1 extrapolate_t 1]);
	tmp2 = repmat(proj(:, end, :), [1 extrapolate_t 1]);
	proj = cat(2, tmp1, proj, tmp2); % [ns, nt+2*extrapolate_t, na]
	nt = nt + 2 * extrapolate_t;
	if offset_t ~= 0
		fail 'offset_t not yet done for extrapolation'
	end
end

switch use_mex
case 1 % loop in mex
	proj = permute(proj, [2 1 3]); % sta -> tsa
	img = back_call('fdk,ts,back', ...
			int32([nx ny nz]), [dx dy dz], ...
			offset_xyz, uint8(mask), ...
			dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
			proj(:,:,ia_list), double(betas(ia_list)), ...
			double(source_zs(ia_list)), ... % DK
			na1_half, ... % DK
			int32(cone_par), int32(w3d), ... % DK
			source_dz_per_rad, single(pitch), ... % DK
			nthread);
	img = permute(img, [2 3 1]); % zxy -> xyz

case 3 % fdk,st (for testing only - slower!)
	warn 'fdk,st does not handle projection view edges completely'
	warn 'fdk,st does not handle helical geometry' % DK 
	img = back_call('fdk,st,back', ...
			int32([nx ny nz]), [dx dy dz], ...
			offset_xyz, uint8(mask), ...
			dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
			proj(:,:,ia_list), double(betas(ia_list)), nthread);

otherwise
	fail 'bug'
end

if scale_dang % final "\der angle" scale:
	warn 'todo: scaling for end slices in mex version'
	%img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
	wz = (nz-1)/2 + offset_xyz(3);
	zc = ([0:nz-1] - wz) * dz;
	for iz=1:nz
		if pitch == 0
			ia_min = 1;
			ia_max = na;
		else
			zmid = zc(iz); % z coordinate of center of this slice
			ia1_middle = imin(abs(zmid - source_zs)); % ''1 because matlab indexing
			ia_min = max(ia1_middle - na1_half, 1);
			ia_max = min(ia1_middle + (na1_half - 1), na);
		end

		orbit_iz = double(((ia_max - ia_min + 1)/ia_skip) * abs(orbit) / (na/ia_skip));
		scale = double(pi / deg2rad(abs(orbit_iz)));
		img(:,:,iz) = double(scale * deg2rad(abs(orbit)) / (na/ia_skip)) * img(:,:,iz);
	end
end

end % cbct_back_mex()


% cbct_back_mat()
% matlab backprojector (slower)
function img = cbct_back_mat(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	source_zs, ...
	source_dz_per_view, pitch, ... % DK
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, scale_dang, extrapolate_t, ...
	cone_par, w3d) % DK

if ~extrapolate_t
	warn 'mat version always extrapolates in t by replication'
	extrapolate_t = true;
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

% DK
na1 = na / orbit * 360;
if abs(round(na1) - na1) > 1e-4
	fail 'bad na/orbit'
end
na1_half = ceil(na1 / 2);
source_dz_per_rad = source_dz_per_view / (deg2rad(abs(orbit)) / (na/ia_skip));

% loop over slices
img = zeros([size(mask) nz]);
sdim = [ns+1 nt]; % trick: extra zeros for zero extrapolation in "s"
proj1 = zeros(sdim);
ticker reset
for iz=1:nz

	if pitch == 0
		ia_min = 1;
		ia_max = na;
	else
	% for helix case, determine relevant view angles by DK (from rebin_helix.m)
		zmid = zc(iz); % z coordinate of center of this slice
		ia1_middle = imin(abs(zmid - source_zs)); % ''1 because matlab indexing
		ia_min = max(ia1_middle - na1_half, 1);
		ia_max = min(ia1_middle + (na1_half - 1), na);
		% todo: for helix case, determine relevant view angles (cf GH version)
	end

	% loop over each projection angle
	img2 = 0;
	for ia=ia_min:ia_skip:ia_max
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc * cos(beta) + yc * sin(beta);
		y_betas = dso - (-xc * sin(beta) + yc * cos(beta));

		% detector indices
		if cone_par % cone_parallel by DK
			mag = dsd ./ (sqrt(dso^2 - x_beta.^2) - (dso - y_betas)); % [np] 
		elseif isinf(dsd) || isinf(dso) % par
			mag = ones(size(y_betas)); % [np]
		else
			mag = dsd ./ y_betas; % [np]
		end

		if cone_par % cone_parallel by DK
			sprime = x_beta; 
		elseif isinf(dfs) ... % flat
			|| isinf(dsd) || isinf(dso) % par
			sprime = mag .* x_beta; 
		elseif (dfs == 0) % arc
			r_loop = x_beta - offset_source;
			sprime = dsd * atan2(r_loop, y_betas);
		end

		% helical by DK
		% same for both arc and flat!
		if cone_par % DK
			source_zs_offset = asin(- x_beta ./ dso) * source_dz_per_rad;
			tprime = mag .* (zc(iz) - (source_zs(ia) + source_zs_offset)); 
		else
			tprime = mag * (zc(iz) - source_zs(ia)); % \tbxyz
		end

		bs = sprime / ds + ws;
		bt = tprime / dt + wt;

		bs(bs < 1 | bs > ns) = ns+1; % trick for zero extrapolation in s
		bt = max(bt, 1); % implicit extrapolation in t (detector row)
		bt = min(bt, nt);

		% bi-linear interpolation:
		is = floor(bs); % left bin
		it = floor(bt);

		is(is == ns+1) = ns; % trick for zero extrapolation in s
		it(it == nt) = nt - 1; % trick for last row extrapolation in t

		wr = bs - is;	% right weight
		wl = 1 - wr;	% left weight
		wu = bt - it;	% upper weight
		wd = 1 - wu;	% lower weight

		proj1(1:ns,:) = proj(:,:,ia); % trick: 1 extra zero at ns+1

		p1 =	wl .* proj1(sub2ind(sdim, is  , it  )) + ...
			wr .* proj1(sub2ind(sdim, is+1, it  ));
		p2 =	wl .* proj1(sub2ind(sdim, is  , it+1)) + ...
			wr .* proj1(sub2ind(sdim, is+1, it+1));
		p0 = wd .* p1 + wu .* p2; % vertical interpolation

		if cone_par % DK
			if w3d
				if pitch 
					kk = w3d * pitch;
					ia_c = max(min(ia + na1_half ...
						* ((tprime >= 0) - (tprime < 0)), na), 1);
				else 
					kk = w3d * abs(zc(iz));
					ia_c = mod(ia + na1_half - 1, na) + 1;
				end
				x_beta_c = - x_beta; 
				y_betas_c = 2 * dso - y_betas; 
				mag_c = dsd ./ (sqrt(dso^2 - x_beta_c.^2) ... 
					- (dso - y_betas_c)); 
				source_zs_offset_c = asin(- x_beta_c ./ dso) ...
						* source_dz_per_rad;
				tprime_c = mag_c .* (zc(iz) ...
					- (source_zs(ia_c) + source_zs_offset_c));
				scale = abs(tprime / dsd).^kk;
				scale_c = abs(tprime_c / dsd).^kk;
			
				%	
				eps = 0.001;
				ww1 = 2 * Min(Max(scale_c ./ (scale_c + scale), eps), 1-eps);
				% no view weighing for rays without their conjugates.
				ww1(ia_c == 1) = 1;
				ww1(ia_c == na) = 1;
				
				p0 = p0 .* ww1;
			end
		elseif isinf(dfs) ... % flat
			|| isinf(dsd) || isinf(dso) % par
			p0 = p0 .* mag.^2; % back-projection weighting for flat
		elseif dfs == 0 % arc
			p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + y_betas.^2);
		end
		
		img2 = img2 + p0;
	end % ia

	% DK
	if scale_dang % final ''\der angle'' scale:
		orbit_iz = ((ia_max - ia_min + 1)/ia_skip) * abs(orbit) / (na/ia_skip);
		scale = pi / deg2rad(abs(orbit_iz));
		img2 = (scale * deg2rad(abs(orbit)) / (na/ia_skip)) * img2;
	end	

	img(:,:,iz) = embed(img2, mask);
end % iz

% DK
%if scale_dang % final "\der angle" scale:
%	img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
%end

end % cbct_back_mat()


% fdk_mex_call()
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


% cbct_back_test()
% compare various versions to check consistency
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

if 0 % check extrapolator for mat version (pointless because always applied!)
	tmp1 = cbct_back(cg.ones, args{2:end}, 'use_mex', 0);
	tmp2 = cbct_back(cg.ones, args{2:end}, 'use_mex', 0, 'extrapolate_t', 1);
	im plc 1 3
	im(1, tmp1)
	im(2, tmp2)
	im(3, tmp2 - tmp1)
return
end

back0 = cbct_back(args{:}, 'use_mex', 0); % mat version as reference
% im(back0), return

if exist('jf_mex') == 3

	if 0 % check extrapolator for mex version
		bmex0 = cbct_back(args{:}, 'use_mex', 1, 'extrapolate_t', 0);
		bmex1 = cbct_back(args{:}, 'use_mex', 1, ...
				'extrapolate_t', ceil(1.3 * cg.nt/2));
		im plc 2 3
		im(1, back0)
		im(2, bmex0)
		im(3, bmex1)
		im(5, bmex0 - back0)
		im(6, bmex1 - back0)
		max_percent_diff(back0, bmex0)
		max_percent_diff(back0, bmex1)
		pr 'nrms(bmex0(:), back0(:))'
		pr 'nrms(bmex1(:), back0(:))'
	return
	end

	for ii = 1:3
		use_mex = ii;
		printm('testing jf_mex with use_mex=%d', use_mex)
		back1{ii} = cbct_back(args{:}, ...
			'extrapolate_t', ceil(1.3 * cg.nt/2), ...
			'use_mex', use_mex, 'back_call', @jf_mex);
	end
	printm('mpd mex1 vs mat: %g%%', max_percent_diff(back1{1}, back0))
	printm('mpd mex1 vs mex2: %g%%', max_percent_diff(back1{1}, back1{2}))
	if f.compare_st % only if edges are zeroed
	printm('mpd mex1 vs mex3: %g%%', max_percent_diff(back1{1}, back1{3}))
	end

	if 1 % examine discrepancies between matlab and mex
		im plc 2 2
		iz = 1:ig.nz;
		im(1, back0(:,:,iz), 'mat'), cbar
		im(2, back1{1}(:,:,iz), 'mex1'), cbar
		err = back1{1}(:,:,iz) - back0(:,:,iz);
		im(3, err, [-1 1]*100), cbar
		im(3, err, [-1 1]*100), cbar
		im(4, jf_mip3(abs(err)), 'mip |diff|'), cbar
%		im(4, jf_mip3(err, 'type', 'mid'), 'diff'), cbar
	end

%	im clf, im_toggle(back0(:,:,iz), back1{1}(:,:,iz))
end

if exist('fdk_mex') == 3
	printm('found fdk_mex:\n %s', which('fdk_mex'))
	% todo: need to consider what to do about extrapolate_t for fdk_mex!
	back2 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @fdk_mex);
	printm('mpd jf_mex vs fdk_mex: %g%%', max_percent_diff(back1{1}, back2))

	if 1
		im plc 1 3
		iz = 1:ig.nz;
		im(1, back1{1}(:,:,iz)), cbar
		im(2, back2(:,:,iz)), cbar
		err = back2(:,:,iz) - back1{1}(:,:,iz);
		im(3, err), cbar
	end
end

end % cbct_back_test()
