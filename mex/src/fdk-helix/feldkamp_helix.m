  function img = feldkamp_helix(cg, ig, proj, varargin)
%|function img = feldkamp_helix(cg, ig, proj, varargin)
%|
%| FBP reconstruction of cone-beam tomography data collected with
%| a circular and ''helical'' source trajectory.
%| See ir_fdk_helix_example.m for example.
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns nt na]	cone-beam projection views (line integrals)
%|
%| options
%|	'window' [npad]		'ramp' (default), or 'hann', or array.
%|				if array, then use samples [-K/2, K/2).
%|	'offset_source'		distance from isocenter to perpendicular ray
%|				[the same units (e.g., mm) as pixel_size etc.]
%|				caution: probably should not be used
%|	'ia_skip' [int]		downsample in angle to save time for tests
%|	'use_mex'	0|1|?	mex backprojector? see cbct_back.m
%|	'nthread'		default: jf('ncore')
%|	'extrapolate_t'	0|?	if 0, no detector row replication;
%|				if >0, repeat top and bottom detector rows
%|				by this many extra rows.
%|				default=0 but recommend using ~ nt/2
%|	'cone_par'		1 if cone-parallel geometry
%|	'w3d'			@tang:05:atd, tang:06:atd
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|
%| References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
%| Notation here follows Fessler tomography book chapter (ask if interested).
%|
%| Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
%| and Jeff Fessler, University of Michigan
%| 2005-06-13 arc detector case contributed by Yingying Zhang
%| 2013-03-20 weighting for arc detector case corrected by Rebecca Malinas,
%| thanks to suggestion by Jinyi Qi.
%| 2014-08-04 changed w1cyl to use proper w1 weighting by default
%|
%| 2014-11-10 helical version of feldkamp.m with view weighting
%| proposed in @tang:05:atd and @tang:06:atd, contributed by Donghwan Kim
%|

if nargin == 1 && streq(cg, 'test')
	run_mfile_local('feldkamp_example')
return
end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.use_mex = []; % defer to cbct_back.m
arg.window = 'ramp'; % 'ramp'
arg.offset_source = 0; % r_off, distance between rotation iso-center
			% and ray from source that is orthogonal to detector.
arg.w1cyl = 1; % use proper w1 weighting that works "perfectly" for cylinders
arg.ia_skip = 1;
arg.extrapolate_t = 0; % no extrapolation by default for historical reasons
arg.nthread = jf('ncore');
arg.cone_par = 0; % DK
arg.w3d = 0; % DK
arg = vararg_pair(arg, varargin);

img = feldkamp_do(proj, ...
	cg, ig, ...
	cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
	cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
	ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
	[ig.offset_x ig.offset_y ig.offset_z], ...
	arg.w1cyl, arg.window, arg.ia_skip, arg.extrapolate_t, ...
	arg.use_mex, arg.nthread, ...
	arg.cone_par, arg.w3d); % DK
end % feldkamp()


% feldkamp_do()
function img = feldkamp_do(proj, ...
	cg, ig, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	w1cyl, window, ia_skip, extrapolate_t, use_mex, nthread, ...
	cone_par, w3d) % DK

% step 1: weight each projection view like in fan-beam case
proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl, cone_par);

% step 2: filter each projection view
[ns nt na] = size(proj);
proj = fdk_filter(proj, window, dsd, dfs, ds);

% step 3: cone-beam backprojection of the filtered views
cpu etic
img = cbct_back_helix(proj, cg, ig, ...
	'offset_source', offset_source, ...
	'ia_skip', ia_skip, ...
	'extrapolate_t', extrapolate_t, ...
	'use_mex', use_mex, ...
	'nthread', nthread, ...
	'cone_par', cone_par, 'w3d', w3d); % DK: todo in this mfile
cpu etoc 'fdk backprojection time:'

end % feldkamp_do()


% feldkamp_weight1()
% step 1: weight the projections as in fan-beam case
function proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, ...
	dsd, dso, dfs, w1cyl, cone_par);
[ns nt na] = size(proj);
ss = ([-(ns-1)/2:(ns-1)/2]' - offset_s) * ds;
tt = ([-(nt-1)/2:(nt-1)/2]' - offset_t) * dt;

[ss tt] = ndgrid(ss, tt);
if cone_par % cone-parallel
	ww1 = dsd ./ sqrt(dsd^2 + tt.^2);
elseif isinf(dfs) % flat
	if w1cyl % weighting that is "exact" for cylindrical-like objects
		ww1 = dso ./ sqrt(dsd^2 + ss.^2 + tt.^2);
	else
		warn 'using old and wrong version of w1 weighting!'
		ww1 = dso ./ sqrt(dsd^2 + ss.^2 + tt.^2) ...
			.* sqrt(1 + (tt/dsd).^2); % original version in book
	end
elseif dfs == 0 % arc
	if w1cyl % weighting that is "exact" for cylindrical-like objects
		ww1 = (dso/dsd) * cos(ss/dsd) ./ sqrt(1 + (tt/dsd).^2);
	else
		warn 'using old and wrong version of w1 weighting!'
		ww1 = (dso/dsd) * cos(ss ./ (dsd * sqrt(1 + (tt/dsd).^2))); % wrong!
	end
else
	error 'other configurations not implemented'
end

for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia) .* ww1;
end
end % feldkamp_weight1()
