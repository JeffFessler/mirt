 function st = ct_geom(type, varargin)
%function st = ct_geom(type, varargin)
%|
%| Create the "CT geometry" structure that describes the sampling
%| characteristics of a cone-beam CT system (axial or helical).
%| (Use sino_geom() for 2D fan-beam or parallel-beam systems.)
%|
%| in
%|	type	'fan'		multi-slice fan-beam - recommended
%|		'par'		(parallel-beam only weakly supported)
%|
%| options for all geometries
%|	'orbit_start'		default: 0
%|	'orbit'			[degrees] default: 180 for parallel / mojette
%|					or 360 for fan (negative for CW)
%|					can be 'short' for fan-beam short scan
%|	'down'			down-sampling factor, for testing
%|					can be [down_s down_t down_a]
%|	'units'			string to print distance units (default: 'mm')
%|
%| options for fan-beam
%|	'ns'			# of horizontal samples
%|	'nt'			# of vertical samples
%|	'na' | 'nbeta'		# of angular samples
%|	'ds'			horizontal sample spacing (default: 1)
%|	'dt'			vertical sample spacing (default: -ds)
%|				or {'dz', dz} to use dz * dsd / dso (usual CT)
%|	'offset_s'		unitless fraction of a channel (default: 0)
%|				(relative to line between two center channels).
%|				use 0.25 or 1.25 for "quarter-detector offset"
%|	'offset_t'		unitless (default: 0)
%|
%| options for partial scans (added by jang-hwan cho)
%|	'nframe'		# of frames to divide orbit into (default: 1)
%|	'frame'			which frame? (default: 1)
%|
%| options for helical (or step-and-shoot)
%|	'pitch'			bed_travel_per_rotation / axial_fov. default: 0.
%|				(unitless. can be negative. usually near 1.0)
%|	'source_z0'		z-location of source for first view. Default 0.
%|					It must have same units as dt and ds.
%|				Use 'center' to center the helix around z=0.
%|
%|	'user_source_zs' [na]	user-specified source z-locations for each view.
%|				usually this is empty (default) in which case
%|				source_zs is computed internally from "pitch"
%|				It must have same units as dt and ds.
%|
%|	fan beam distances:
%|	'dsd' | 'dis_src_det'	default: inf (parallel beam)
%|	'dso' | 'dis_src_iso'	default: inf (parallel beam)
%|	'dod' | 'dis_iso_det'	default: 0
%|	'dfs' | 'dis_foc_src'	default: 0 (3rd generation CT arc),
%|					use 'inf' for flat detector
%|
%| out
%|	st	(struct)	initialized structure
%|
%| methods
%|	st.shape(sino)		reshape sinograms that are columns into 3d array
%|	st.s			s sample locations
%|	st.t			t sample locations
%|	st.gamma		[nb] gamma sample values [radians]
%|	st.gamma_s		gamma values given s values
%|	st.gamma_max		half of fan angle [radians] (if offset_s=0)
%|	st.gamma_max_abs	half of fan angle [radians] (recommended)
%|	st.ws			(ns-1)/2 + st.offset_s
%|	st.wt			(nt-1)/2 + st.offset_t
%|	st.ad			[na] source angles in degrees
%|	st.ar			[na] source angles in radians
%|	st.dim			dimensions: [st.ns st.nt st.na]
%|	st.downsample(down)	reduce sampling by integer factor
%|	st.ones			ones(ns,nt,na, 'single')
%|	st.zeros		zeros(ns,nt,na, 'single')
%|	st.rad			[ns] radial distance of each ray from origin
%|	st.rmax			max radius within FOV
%|	st.footprint_size(ig)	max footprint width in 's'
%|	st.cone_angle		(half) cone angle on axis (s=0): +/- angle
%|	st.zfov			axial FOV
%|	st.source_zs		[na] z-locations of source for each view
%|	st.shape(sino(:))	reshape to [ns,nt,na,?]
%|	st.unitv(is,it,ia)	unit 'vector' with one nonzero element
%|	st.xds st.yds		x,y center coordinates of middle detector row
%|	st.xd_s(s) st.yd_s(s)	x,y center coordinates of middle detector row given s
%|	st.plot([ig])		show geometry
%|
%|	trick: you can make orbit=0 and orbit_start = column vector (length na)
%|	if you need nonuniformly spaced projection view angles.
%|
%| Copyright 2006-1-18, Jeff Fessler and Jang-Hwan Cho, University of Michigan
%|
%| 2009-12-04 modified source_zs definition to use source_z0, eliminate na/2

if nargin == 1 && streq(type, 'test'), ct_geom_test, return, end
if nargin < 1, ir_usage, end

if streq(type, 'ge1') % special case: GE fan-beam
	st = ct_geom_ge1(type, varargin{:});
return
end

if streq(type, 'ge2') % special case: GE axial or helical
	st = ct_geom_ge2(type, varargin{:});
return
end

if streq(type, 'hd1') % special case: GE HD (UM only)
	st = ct_geom_hd1(type, varargin{:});
return
end

if streq(type, 'revo1-cone') % special case: GE Revo1 (UM only)
	tmp = ir_fan_geom_revo1(type);
	st = ct_geom('fan', tmp{:}, 'pitch', 0, 'source_z0', 0, varargin{:});
return
end

% defaults
st.type = type;
st.ns = [];
st.nt = [];
st.na = [];
st.down = 1;
st.nframe = 1; % entire scan as 1 "frame"
st.frame = 1;
st.orbit_start = 0;
st.pitch = 0; % default for axial
st.source_z0 = 0; % default for axial
st.units = 'mm';
st.user_source_zs = [];

if streq(type, 'fan')
	st = ct_geom_fan(st, varargin{:});
%elseif streq(type, 'par')
%	st = ct_geom_par(st, varargin{:});
%elseif streq(type, 'moj')
%	st = ct_geom_moj(st, varargin{:});
else
	fail('unknown sinotype %s', type)
end

if isempty(st.na), st.na = 2 * floor(st.ns * pi/2 / 2); end

ct_geom_checks(st)

meth = { ...
	's', @ct_geom_s, '()'; ...
	't', @ct_geom_t, '()'; ...
	'ws', @ct_geom_ws, '()'; ...
	'wt', @ct_geom_wt, '()'; ...
	'ad', @ct_geom_ad, '()'; ...
	'ar', @ct_geom_ar, '()'; ...
	'gamma', @ct_geom_gamma, '()'; ...
	'gamma_s', @ct_geom_gamma_s, '(s)'; ...
	'gamma_max', @ct_geom_gamma_max, '()'; ...
	'gamma_max_abs', @ct_geom_gamma_max_abs, '()'; ...
	'cone_angle', @ct_geom_cone_angle, '()'; ...
	'zfov', @ct_geom_zfov, '()'; ...
	'source_dz_per_view', @ct_geom_source_dz_per_view, '()'; ...
	'source_zs', @ct_geom_source_zs, '() -> [na]'; ...
	'orbit_short', @ct_geom_orbit_short, '()'; ...
	'xd_s', @ct_geom_xd_s, '(s)'; ...
	'yd_s', @ct_geom_yd_s, '(s)'; ...
	'xds', @ct_geom_xds, '()'; ...
	'yds', @ct_geom_yds, '()'; ...
	'downsample', @ct_geom_downsample, '()'; ...
	'dim', @ct_geom_dim, '()'; ...
	'ones', @ct_geom_ones, '()'; ...
	'rad', @ct_geom_rad, '()'; ...
	'rmax', @ct_geom_rmax, '()'; ...
	'footprint_size', @ct_geom_footprint_size, '(ig)'; ...
	'unitv', @ct_geom_unitv, '() | (is,it,ia)'; ...
	'zeros', @ct_geom_zeros, '()'; ...
	'shape', @ct_geom_shape, '()'; ...
	'plot', @ct_geom_plot, '() | (ig)';
	'plot3', @ct_geom_plot3, '() | (ig)';
	};

st = strum(st, meth);

if streq(st.source_z0, 'center')
	st.source_z0 = -st.na/2 * st.source_dz_per_view;
end

if any(st.down ~= 1)
	down = st.down; st.down = 1; % trick
	st = st.downsample(down);
end

if streq(type, 'fan') && streq(st.orbit, 'short')
	st.orbit = st.orbit_short / st.nframe; % jc
	st.orbit_start = st.orbit_start + (st.frame - 1) * st.orbit; % jc/jf
	if st.frame ~= 1, warn 'untested and possibly incorrect', end % jf
end


% ct_geom_orbit_short()
function os = ct_geom_orbit_short(st)
os = 180 + 2 * rad2deg(st.gamma_max);


% ct_geom_dim()
function dim = ct_geom_dim(st)
dim = [st.ns st.nt st.na];
if isempty(st.ns) || isempty(st.nt) || isempty(st.na)
	fail 'dim requested without ns,nt,na'
end


% ct_geom_ones()
% sinogram of all ones
function out = ct_geom_ones(st)
out = ones(st.dim, 'single');


% ct_geom_zeros()
% sinogram of all zeros
function out = ct_geom_zeros(st)
out = zeros(st.dim, 'single');


% ct_geom_unitv()
% sinogram where just one ray value is unity (the rest are zero)
function out = ct_geom_unitv(st, is, it, ia)
out = st.zeros;
if ~isvar('is') || isempty(is)
	is = floor(st.ns/2 + 1);
	it = floor(st.nt/2 + 1);
	ia = 1;
end
out(is,it,ia) = 1;


% ct_geom_rad()
% radial distance from origin for each ray
function rad = ct_geom_rad(st)
switch st.type
case 'par'
	rad = st.s
case 'fan'
	if isinf(st.dso) % parallel
		rad = st.s;
	elseif st.dfs == 0 % arc
		gamma = ct_geom_gamma(st);
		rad = st.dso * sin(gamma);
	elseif isinf(st.dfs) % flat
		gamma = ct_geom_gamma(st);
		rad = st.dso * sin(gamma);
	else
		fail 'unknown case'
	end
otherwise
	fail 'unknown type'
end


% ct_geom_rmax()
% max radius within fov
function rmax = ct_geom_rmax(st)
smax = max(abs(st.s));
if streq(st.type, 'fan')
	if isinf(st.dso) % parallel
		rmax = smax;
	elseif st.dfs == 0 % arc
		gamma_max_abs = ct_geom_gamma_max_abs(st);
		rmax = st.dso * sin(gamma_max_abs);
%		rmax2 = st.dso * sin(smax / st.dsd); % todo: cut
%		jf_equal(rmax, rmax2)
	elseif isinf(st.dfs) % flat
		gamma_max_abs = ct_geom_gamma_max_abs(st);
		rmax = st.dso * sin(gamma_max_abs);
%		rmax2 = st.dso * sin(atan(smax / st.dsd)); % todo: cut
%		jf_equal(rmax, rmax2)
	else
		fail 'unknown case'
	end
end


% ct_geom_ws()
% 'middle' sample position
function ws = ct_geom_ws(st)
ws = (st.ns-1)/2 + st.offset_s;


% ct_geom_wt()
% 'middle' sample position
function wt = ct_geom_wt(st)
wt = (st.nt-1)/2 + st.offset_t;


% ct_geom_s()
% sample locations ('radial')
function s = ct_geom_s(st, varargin)
s = st.ds * ([0:st.ns-1]' - st.ws);
if length(varargin)
	s = s(varargin{:});
end


% ct_geom_t()
% sample locations ('radial')
function t = ct_geom_t(st, varargin)
t = st.dt * ([0:st.nt-1]' - st.wt);
if length(varargin)
	t = t(varargin{:});
end


% ct_geom_ad()
% angular sample locations (degrees)
function ang = ct_geom_ad(st, varargin)
ang = [0:st.na-1]'/st.na * st.orbit + st.orbit_start;
ang = ang(varargin{:});

% ct_geom_ar()
% angular sample locations (radians)
function ang = ct_geom_ar(st, varargin)
ang = deg2rad(ct_geom_ad(st));
ang = ang(varargin{:});


% ct_geom_shape()
% reshape into sinogram array
function sino = ct_geom_shape(st, sino)
sino = reshapee(sino, st.ns, st.nt, st.na, []);


% ct_geom_downsample()
% down-sample (for testing)
function st = ct_geom_downsample(si, down)
st = si;
st.down = st.down .* down;

switch numel(down)
case 1
	down_s = down;
	down_t = down;
	down_a = down;
case 3
	down_s = down(1);
	down_t = down(2);
	down_a = down(3);
otherwise
	fail('bad down %g', down)
end
clear down

st.ns = 4 * ceil(st.ns / down_s / 4); % multiple of 4 for simd
st.nt = 2 * ceil(st.nt / down_t / 2); % keep it even

user_specified = false; % did user specify zs or orbit_start vector?

if ~isempty(st.user_source_zs)
	st.user_source_zs = st.user_source_zs(1:down_a:st.na);
	user_specified = true;
end

if numel(st.orbit_start) > 1
	st.orbit_start = st.orbit_start(1:down_a:si.na);
	user_specified = true;
end

if streq(st.type, 'fan')
	st.ds = st.ds * down_s;
	st.dt = st.dt * down_t;
else
	fail('unknown sinotype "%s"', type)
end

if user_specified % if user-specified then just decimate
	st.na = length([1:down_a:st.na]);

elseif all(diff(si.source_zs) == 0) % axial
	st.na = max(1, round(si.na / down_a)); % at least one view

else % helical described by pitch and source_z0
	nturn = si.orbit / 360;
	na1 = si.na / nturn; % # of views in 1 turn, originally

	% if it is helical with equal view spacing, then adjust orbit slightly,
	% to preserve integer number of views per turn (if applicable)
	tol = 1e-4;
	if abs(round(na1) - na1) < tol % integer views per turn
		na2 = round(na1 / down_a); % new # of views in 1 turn (integer)
		tmp = nturn * na2; % new na value, roughly na/down
		if abs(round(tmp) - tmp) < tol
			st.na = round(tmp);
		else % adjust orbit
			st.na = floor(tmp);
			st.orbit = 360 / na2 * st.na;
		end
	else
		st.na = round(si.na / down_a);
	end
end


% ct_geom_cone_angle()
function cone_angle = ct_geom_cone_angle(st)
if isinf(st.dso) || isinf(st.dsd) % parallel beam
	cone_angle = 0;
else % cone
	cone_angle = atan((st.nt * st.dt)/2 / st.dsd);
end


% ct_geom_zfov()
% axial FOV, considering magnification factor, as collimated at iso
function zfov = ct_geom_zfov(st)
if isinf(st.dso) || isinf(st.dsd) % parallel beam
	zfov = st.nt * st.dt;
else % cone
	zfov = st.dso / st.dsd * st.nt * st.dt;
end


% ct_geom_source_dz_per_view()
% source axial travel for each view
function out = ct_geom_source_dz_per_view(st)
if ~isempty(st.user_source_zs), fail 'undefined', end
if st.na == 1 % for tomosynthesis
	out = 0;
	return
end
if ~st.pitch
	out = 0;
	return
end
if length(st.orbit) ~= 1 || st.orbit == 0
	fail 'dz_per_view undefined for user angles'
end
na_per_360 = st.na * (360 / st.orbit); % # views per turn
out = st.pitch * st.zfov / na_per_360;


% ct_geom_checks
% check vectors user_source_zs and orbit_start
function ct_geom_checks(st)
if ~isempty(st.user_source_zs) && (st.pitch ~= 0)
	fail('only one of "pitch" (recommended) or user_source_zs can be given')
end
if ~isempty(st.user_source_zs) && length(st.user_source_zs) ~= st.na
	fail('user_source_zs size mismatch')
end
if numel(st.orbit_start) > 1 && numel(st.orbit_start) ~= st.na
	fail('orbit_start vector length')
end
if numel(st.orbit_start) > 1 && (st.pitch ~= 0) && isempty(st.user_source_zs)
	fail('helical with vector orbit_start: must specify user_source_zs')
end
if round(st.ns) ~= st.ns, fail('ns must be integer'), end
if round(st.nt) ~= st.nt, fail('nt must be integer'), end
if round(st.na) ~= st.na, fail('na must be integer'), end


% ct_geom_source_zs()
% source z locations
function source_zs = ct_geom_source_zs(st, varargin)

if ~isempty(st.user_source_zs) % user-provided
	source_zs = st.user_source_zs;
	source_zs = source_zs(varargin{:});
else
	source_dz = st.source_dz_per_view;
	source_zs = st.source_z0 + [0:st.na-1]' * source_dz;
	source_zs = source_zs(varargin{:});
end


% ct_geom_gamma_s()
% gamma values from s values
function gamma = ct_geom_gamma_s(st, ss)
switch st.dfs
case 0
	gamma = ss / st.dsd; % 3rd gen: equiangular arc
case inf
	gamma = atan(ss / st.dsd); % flat
otherwise
	fail 'dfs not done'
end


% ct_geom_gamma()
% gamma sample values
function gamma = ct_geom_gamma(st, varargin)
gamma = ct_geom_gamma_s(st, st.s);
gamma = gamma(varargin{:});


% ct_geom_gamma_max()
function gamma_max = ct_geom_gamma_max(st)
gamma_max = max(st.gamma);


% ct_geom_gamma_max_abs()
function gamma_max = ct_geom_gamma_max_abs(st)
gamma_max = max(abs(st.gamma));


% ct_geom_xd_s()
% x positions along detector given s values (for beta=0)
function xds = ct_geom_xd_s(st, ss)
if nargin < 2
	ss = st.s;
end
switch st.type
case 'par'
	xds = ss;
case 'fan'
	switch st.dfs
	case 0 % arc
		xds = st.dsd * sin(st.gamma_s(ss));
	case inf % flat
		xds = ss;
	otherwise
		fail 'not done'
	end
otherwise
	fail 'bug'
end


% ct_geom_yd_s()
% y positions along detector given s values (for beta=0)
function yds = ct_geom_yd_s(st, ss)
if nargin < 2
	ss = st.s;
end
switch st.type
case 'par'
	yds = zeros(size(ss), class(ss));
case 'fan'
	switch st.dfs
	case 0 % arc
		yds = st.dso - st.dsd * cos(st.gamma_s(ss));
	case inf % flat
		yds = -st.dod * ones(size(ss), class(ss));
	otherwise
		fail 'not done'
	end
otherwise
	fail 'bug'
end


% ct_geom_xds()
% x center positions of detectors (for beta=0)
function xds = ct_geom_xds(st, varargin)
xds = ct_geom_xd_s(st, st.s);
xds = xds(varargin{:});


% ct_geom_yds()
% y center positions of detectors (for beta=0)
function yds = ct_geom_yds(st, varargin)
yds = ct_geom_yd_s(st, st.s);
yds = yds(varargin{:});


% ct_geom_footprint_size()
% biggest possible footprint over image
function footprint_size = ct_geom_footprint_size(st, ig, varargin)
di = sqrt(ig.dx^2 + ig.dy^2);
smax = max(abs(st.s));
rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2;
dso = st.dso;
dsd = st.dsd;
if rfov > 0.99 * dso, fail 'bad dso', end

if isinf(st.dso) % parallel
	footprint_size = di;
elseif st.dfs == 0 % arc3
	footprint_size = di * dsd / (dso - rfov);
elseif isinf(st.dfs) % flat
	footprint_size = di * sqrt(dsd^2 + smax^2) / (dso - rfov);
else
	fail 'unknown case'
end
footprint_size = footprint_size / st.ds; % unitless 


%
% ct_geom_fan()
%
function st = ct_geom_fan(st, varargin);

% defaults
st.orbit = 360; % [degrees]
st.ds		= 1;
st.dt		= [];
st.offset_s	= 0;
st.offset_t	= 0;

st.dsd = [];	% dis_src_det
st.dso = [];	% dis_src_iso
st.dod = [];	% dis_iso_det
st.dfs = 0;	% dis_foc_src (3rd gen CT)

subs = { ...
	'src_det_dis', 'dsd';
	'dis_src_det', 'dsd';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	'nbeta', 'na';
	'source_zs', 'user_source_zs';
	};
st = vararg_pair(st, varargin, 'subs', subs);

% work out distances
if (~isempty(st.dsd) && isinf(st.dsd)) ...
|| (~isempty(st.dso) && isinf(st.dso)) % handle parallel-beam case gracefully
	st.dsd = inf; st.dso = inf; st.dod = 1;
end
if isempty(st.dsd) + isempty(st.dso) + isempty(st.dod) > 1
	fail 'must provide at least two of dsd, dso, dod'
end
if isempty(st.dsd), st.dsd = st.dso + st.dod; end
if isempty(st.dso), st.dso = st.dsd - st.dod; end
if isempty(st.dod), st.dod = st.dsd - st.dso; end
% caution: we do the following check here, but if later the user changes
% any of the three quantitities alone then we may lose consistency.
% Gtomosyn.m is an example where 'dsd' and 'dso' are changed on the fly.
if st.dso + st.dod ~= st.dsd
	fail 'bad fan-beam distances'
end

if isempty(st.dt), st.dt = -st.ds; end
if iscell(st.dt)
	if length(st.dt) == 2 && streq(st.dt{1}, 'dz')
		dz = st.dt{2};
		st.dt = dz * st.dsd / st.dso;
		printm('dt = dz * dsd / dso = %g * %g / %g = %g', ...
			dz, st.dsd, st.dso, st.dt)
	else
		fail 'bad dt cell'
	end
end


% ct_geom_plot2()
% picture of 2D source position / detector geometry
function ct_geom_plot2(st, ig, varargin)
arg.tomosyn = false; % set to 1 for offset_x plotting trick
arg = vararg_pair(arg, varargin);
if arg.tomosyn, fail 'not done', end
if ~streq(st.type, 'fan'), fail 'only fan done', end
t = linspace(0,2*pi,1001);

if isinf(st.dsd) % parallel beam
	rfov = max(abs(st.s));
	plot(	0, 0, '.', ...
		rfov * cos(t), rfov * sin(t), 'm:', ...
		st.s, -rfov, 'yo')

else % fan beam

	beta0 = deg2rad(st.orbit_start);
	rot = [cos(beta0) sin(beta0); -sin(beta0) cos(beta0)];
	p0 = st.dso * [-sin(beta0); cos(beta0)]; % rot' * [0; st.dso];
	pd = rot' * [st.xds'; st.yds']; % rotated arc detector points
	rfov = st.dso * sin(max(abs(st.gamma)));

	plot(	p0(1), p0(2), 'ro', ... % source
		st.dso * cos(t), st.dso * sin(t), 'c--', ... % source circle
		pd(1,:), pd(2,:), 'yo') % detector elements

	hold on
	if isvar('ig') && ~isempty(ig)
		xmin = min(ig.x); xmax = max(ig.x);
		ymin = min(ig.y); ymax = max(ig.y);
		im(ig.x, ig.y, sum(ig.mask,3))

		plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	end

	plot(	0, 0, '.', ...
		[pd(1,1) p0(1) pd(1,end)], [pd(2,1) p0(2) pd(2,end)], 'r-', ...
		rfov * cos(t), rfov * sin(t), 'm:')

	hold off
	axis tight
	tmp = truncate_precision(1.05*max(axis), 2, 'type', 'ceil');
	axis([-1 1 -1 1]*tmp)
%	axis equal
%	axis square

end


if isvar('ig') && ~isempty(ig)
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	hold off
end
titlef('rfov = %g', rfov)
axis square


% ct_geom_plot3()
% picture of 3D helical CT geometry
function dummy = ct_geom_plot3(st, ig)
dummy = [];
if ~streq(st.type, 'fan'), fail 'only fan done', end

if isinf(st.dso) || isinf(st.dsd)
	warn 'parallel-beam plot not done'
return
end

t1 = -st.dso * sin(st.ar); % source x pos
t2 = st.dso * cos(st.ar); % source y pos
t3 = st.source_zs; % source z pos
plot3(t1, t2, t3, 'y.') % source trajectory
view(-110, 22)
hold on
for ia = [1 1+floor(st.na/2) st.na]
	text(t1(ia), t2(ia), t3(ia), sprintf('%d', ia))
end

% axes
r100 = @(x) 100 * ceil(max(abs(x)) / 100);
r10 = @(x) 10 * ceil(max(abs(x)) / 10);
plot3([-1 1] * r100(t1), [0 0], [0 0])
text(r100(t1), 0, 'x')
plot3([0 0], [-1 1] * r100(t2), [0 0])
text(0, r100(t2), 'y')
plot3([0 0], [0 0], [-1 1] * r10(t3))
text(0, 0, r10(t3), 'z')

grid on
xlabelf('x (%s)', st.units)
ylabelf('y (%s)', st.units)
zlabelf('z (%s)', st.units)
if ~st.pitch % meaningful ticks for axial case
	tmp = [st.zfov/2 ig.dz*ig.nz/2 st.nt*st.dt/2];
	ztick(unique(sort([-tmp 0 tmp])))
end

if isvar('ig') && ~isempty(ig) && isfield(ig.meth,'z') % image box
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	zmin = min(ig.z); zmax = max(ig.z);
	plot3(	[xmin xmax xmax xmin xmin ], ...
		[ymin ymin ymax ymax ymin ], ...
		[zmin zmin zmin zmin zmin ], ...
		'g-')
	plot3(	[xmin xmax xmax xmin xmin], ...
		[ymin ymin ymax ymax ymin], ...
		[zmax zmax zmax zmax zmax], ...
		'g-')
	plot3(xmin*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymin*[1 1], [zmin zmax], 'g-')
	plot3(xmin*[1 1], ymax*[1 1], [zmin zmax], 'g-')
	plot3(xmax*[1 1], ymax*[1 1], [zmin zmax], 'g-')
	plot3(xmin*ones(ig.nz,1), ymin*ones(ig.nz,1), ig.z, 'g.')
end

% detector array for source at each position
detcolor = {'r.', 'c.', 'm.'};
switch st.na
case 1
	ia_list = 1;
case 2
	ia_list = [1 st.na];
otherwise 
	ia_list = [1 1+floor(st.na/2) st.na];
end
for ia = ia_list
	src = [t1(ia) t2(ia) t3(ia)];
	unit_vec = [-src(1) -src(2) 0] / sqrt(src(1)^2+src(2)^2);
	% diametrically opposite point of the source on the detector
	det_cen_loc = src + st.dsd * unit_vec;
	plot3(det_cen_loc(1), det_cen_loc(2), det_cen_loc(3), 'ys')

	rot = st.ar(ia);
	rot = [cos(rot) -sin(rot); sin(rot) cos(rot)];
	pd = rot * [st.xds'; st.yds'];

	for it=1:st.nt
		plot3(pd(1,:), pd(2,:), src(3) + st.t(it)*ones(1,st.ns), ...
			detcolor{find(ia == ia_list)})
	end

	if 1 % line from source to middle of first,last detector row
		plot3([src(1) det_cen_loc(1)], [src(2) det_cen_loc(2)], ...
			[src(3) st.source_zs(ia)+st.t(1)], 'm-')
		plot3([src(1) det_cen_loc(1)], [src(2) det_cen_loc(2)], ...
			[src(3) st.source_zs(ia)+st.t(end)], 'm-')
	end
end
hold off


% ct_geom_plot()
function out = ct_geom_plot(st, varargin)
%im clf
if all(st.source_zs == 0) % 2D or axial case
	ct_geom_plot2(st, varargin{:});
else
	ct_geom_plot3(st, varargin{:});
end
if nargout, out = []; end


% ct_geom_ge1()
% 'lightspeed';
% these numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
function st = ct_geom_ge1(type, varargin)

tmp.orbit = 360;
tmp = vararg_pair(tmp, varargin, 'allow_new', true);
if streq(tmp.orbit, 'short')
	na = 642; % trick: reduce na for short scans
	for ii=1:length(varargin)
		if (streq(varargin{ii}, 'orbit'))
			varargin{ii+1} = na/984*360;
			break
		end
	end
else
	na = 984;
end

st = ct_geom('fan', ...
	'ns', 888, ... % detector channels
	'nt', 1, ... % detector rows
	'na', na, ... % angular samples
	'orbit', 360, ...
	'offset_s', 1.25, ... % quarter-detector offset
	'dsd', 949.075, ...
	'dso', 541, ...
	'dfs', 0, ... % arc
	'ds', 1.0239, ... % detector pitch
	'dt', 1.0964, ... % detector row spacing for 0.625mm slices, 2009-12-06
	varargin{:});
%	'dod', 408.075, ...
% 'strip_width', [], ...


% ct_geom_ge2()
% helical CT
function st = ct_geom_ge2(type, varargin)
st = ct_geom_ge1(type, ...
	'nt', 64, ... % 64-slice
	'dt', 949.0750 / 541 * 0.625, ... % about 1.0964
	varargin{:});


% ct_geom_test()
function ct_geom_test
% axial cone-beam test
cg = ct_geom('fan', 'ns', 888, 'nt', 64, 'na', 984, ...
	'offset_s', 1.25, ...
	'dsd', 949, 'dod', 408);
cg.ad(2);
cg.downsample(2);
cg.ws;
cg.nframe; cg.frame; cg.gamma;
cg.gamma_max; cg.gamma_max_abs;
cg.s;
cg.s(cg.ns/2+1);
jf_equal(cg.gamma_s(cg.s(1:4)), cg.gamma(1:4))
cg.rmax;
cg.zfov;
cg.cone_angle;
if im
	cg.plot;
prompt
end

% test user_source_zs
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'source_zs', []);
cg.source_zs(1:2:7);

% helical cone-beam test
cg = ct_geom('fan', 'dsd', 949, 'dod', 408, ...
	'ns', 888, 'ds', 1.0239, 'offset_s', 1.25, ...
	'nt', 8, 'dt', 1.0964, ...
	'na', 2*984, 'orbit', 2*360, ...
	'source_z0', 'center', 'pitch', 1.0);

cg.source_zs(1:2:7);
cg.downsample(2);
cg.s(cg.ns/2+1);

if im
	ig3 = image_geom('nx', 512, 'fov', 500, ...
		'nz', 16, 'dz', 0.625, 'down', 2);
	cg.plot(ig3);
end

if 1 % footprint
	cg = ct_geom('ge2');
	ig = image_geom('nx', 512, 'fov', 500);
	cg.footprint_size(ig);
end
