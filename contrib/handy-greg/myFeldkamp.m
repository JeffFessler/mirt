  function img = myFeldkamp(cg, ig, proj, varargin)
%|function img = feldkamp(cg, ig, proj, varargin)
%|
%| FBP reconstruction of cone-beam tomography data collected with
%| a circular source trajectory.
%| See feldkamp_example.m for example.
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
%|	'use_mex'	0|1	backprojector: 0 for matlab, 1 for mex (default)
%|	'nthread' 		default: jf('ncore')
%|	'w1cyl'		0|1	1 use w1 weighting for cylinders (default: 0)
%|				0 use w1 weighting in book before 2010-1-21
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|
%| References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
%| Notation here follows Fessler tomography book chapter (ask if interested).
%|
%| Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
%| and Jeff Fessler, University of Michigan
%| arc detector case contributed by Yingying Zhang 2005-6-13

if nargin == 1 && streq(cg, 'test')
	run_mfile_local('feldkamp_example')
return
end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.use_mex = 1;
arg.window = 'ramp';
arg.offset_source = 0; % r_off, distance between rotation iso-center
			% and ray from source that is orthogonal to detector.
arg.w1cyl = 0; % use w1 weighting that works for cylinders?
arg.ia_skip = 1;
arg.nthread = jf('ncore');
arg = vararg_pair(arg, varargin);

%if cg.pitch ~= 0 || any(cg.source_zs ~= 0)
%	fail('sorry, helical CT unsupported')
%end

img = feldkamp_do(proj, ...
	cg, ig, ...
	cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
	cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
	ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
	[ig.offset_x ig.offset_y ig.offset_z], ...
	arg.w1cyl, arg.window, arg.ia_skip, arg.use_mex, arg.nthread);
end % feldkamp()

%
% feldkamp_do()
%
function img = feldkamp_do(proj, ...
	cg, ig, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	w1cyl, window, ia_skip, use_mex, nthread)

% step 1: weight each projection view like in fan-beam case
proj = feldkamp_weight1(proj, cg.pitch, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl);

% step 2: filter each projection view
[ns nt na] = size(proj);

proj = fdk_filter(proj, window, dsd, dfs, ds);
% step 3: cone-beam backprojection of the filtered views
cpu etic
img = my_cbct_back(proj, cg, ig, ...
	'offset_source', offset_source, ...
	'ia_skip', ia_skip, ...
	'use_mex', 0, ...
	'nthread', nthread);
cpu etoc 'fdk backprojection time:'

end % feldkamp_do()


%
% feldkamp_weight1()
% step 1: weight the projections as in fan-beam case
%
function proj = feldkamp_weight1(proj, pitch, ds, dt, offset_s, offset_t, ...
	dsd, dso, dfs, w1cyl);
[ns nt na] = size(proj);
ss = ([-(ns-1)/2:(ns-1)/2]' - offset_s) * ds;
tt = ([-(nt-1)/2:(nt-1)/2]' - offset_t) * dt;

[ss tt] = ndgrid(ss, tt);
if pitch ~= 0
    ww1 = dso ./ sqrt(dsd^2 + ss.^2 + tt.^2);
elseif isinf(dfs) % flat
	ww1 = dso ./ sqrt(dsd^2 + ss.^2 + tt.^2);
	if ~w1cyl % Jinyi Qi suggested not to do this for cylinder objects
		ww1 = ww1 .* sqrt(1 + (tt/dsd).^2); % original version in book
	end
elseif dfs == 0 % arc
	ww1 = (dso/dsd) * cos(ss ./ (dsd * sqrt(1 + (tt/dsd).^2)));
	if w1cyl
		ww1 = ww1 ./ sqrt(1 + (tt/dsd).^2); % todo: new version
	end
else
	error 'other configurations not implemented'
end

for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia) .* ww1;
end
end % feldkamp_weight1()
