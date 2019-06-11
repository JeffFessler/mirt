 function st = sino_geom(type, varargin)
%function st = sino_geom(type, varargin)
%|
%| Create the "sinogram geometry" structure that describes the sampling
%| characteristics of a given sinogram for a 2D parallel or fan-beam system.
%| Using this structure facilitates "object oriented" code.
%| (Use ct_geom() instead for 3D axial or helical cone-beam CT.)
%|
%| in
%|	type	'fan' (fan-beam) | 'par' (parallel-beam) | 'moj' (mojette)
%|
%| options for all geometries
%|	'orbit_start'		default: 0
%|	'orbit'			[degrees] default: 180 for parallel / mojette
%|					or 360 for fan
%|					can be 'short' for fan-beam short scan
%|	'down'			down-sampling factor, for testing
%|	'units'			user-specified string, e.g. 'cm'. default: ''
%|
%| options for parallel-beam
%|	'nb' | 'nr'			# radial samples
%|	'na' | 'nphi'			# angular samples
%|	'dr' | 'ray_spacing'		(default: 1)
%|	'offset_r' | 'channel_offset'	unitless (default: 0)
%|	'strip_width'			'd' or 'dr' to equal dr, (default: [])
%|	('d' is highly recommended; default is 0 for historical reasons only)
%|
%| options for mojette are same as parallel except:
%|	'dx' instead of 'dr'
%|
%| options for fan-beam
%|	'nb' | 'ns'			# "radial" samples (along detector)
%|	'na' | 'nbeta'			# angular samples
%|	'ds' | 'ray_spacing'		(default: 1)
%|	'offset_s' | 'channel_offset'	unitless (default: 0)
%|			(relative to centerline between two central channels).
%|			use 0.25 or 1.25 for "quarter-detector offset"
%|	'source_offset'		same units as dx, ds, etc., e.g., [mm]
%|				use with caution!
%|	'strip_width'			'd' or 'ds' to equal ds, (default: [])
%|
%|	fan beam distances:
%|	'dsd' | 'Dsd' | 'dis_src_det'	default: inf (parallel beam)
%|	'dso' | 'Dso' | 'dis_src_iso'	default: inf (parallel beam)
%|	'dod' |	'Dod' | 'dis_iso_det'	default: 0
%|	'dfs' | 'Dfs' | 'dis_foc_src'	default: 0 (3rd generation CT arc),
%|					use 'inf' for flat detector
%|
%| out
%|	st	(struct)	initialized structure
%|
%| methods
%|	st.shape(sino)		reshape sinograms that are columns into array
%|	st.d			ds or dr
%|	st.s			[nb] s sample locations
%|	st.gamma		[nb] gamma sample values [radians]
%|	st.gamma_max		half of fan angle [radians]
%|	st.w			(nb-1)/2 + st.offset
%|	st.ad			source angles [degrees]
%|	st.ar			source angles [radians]
%|	st.dim			dimensions: [st.nb st.na]
%|	st.downsample(down)	reduce sampling by integer factor
%|	st.offset		offset_s or offset_r
%|	st.ones			ones(nb,na)
%|	st.zeros		zeros(nb,na)
%|	st.rfov			radial fov
%|	st.xds			[nb 1] center of detector elements (beta=0)
%|	st.yds			[nb 1] ""
%|	st.shape(sino(:))	reshape to [nb na ?]
%|	st.unitv(ib,ia)		unit 'vector' with single nonzero element
%|	st.taufun(x,y)		projected s/ds for each (x,y) pair
%|	st.plot			plot system geometry (for fan)
%|
%| Copyright 2006-1-18, Jeff Fessler, University of Michigan

if nargin == 1 && streq(type, 'test'), sino_geom_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if streq(type, 'ge1') % special case
	st = sino_geom_ge1(type, varargin{:});
return
end

if streq(type, 'hd1') % special case
	st = sino_geom_hd1(type, varargin{:});
return
end

if streq(type, 'revo1-fan') % special case
	tmp = ir_fan_geom_revo1(type);
	st = sino_geom('fan', tmp{:}, varargin{:});
return
end

% defaults
st.type = type;
st.nb = [];
st.na = [];
st.down = 1;
st.units = '';
st.orbit_start = 0;
st.source_offset = 0;

meth = { ...
	'd', @sino_geom_d, '()'; ...
	's', @sino_geom_s, '()'; ...
	'w', @sino_geom_w, '()'; ...
	'ad', @sino_geom_ad, '()'; ...
	'ar', @sino_geom_ar, '()'; ...
	'downsample', @sino_geom_downsample, '()'; ...
	'offset', @sino_geom_offset, '()'; ...
	'dim', @sino_geom_dim, '()'; ...
	'ones', @sino_geom_ones, '()'; ...
	'plot', @sino_geom_plot, '() | (ig)'; ...
	'taufun', @sino_geom_taufun, '(x,y) -> [numel(x) na]'; ...
	'unitv', @sino_geom_unitv, '(ib,ia)'; ...
	'zeros', @sino_geom_zeros, '()'; ...
	'shape', @sino_geom_shape, '()'; ...
	'xds', @sino_geom_xds, '()'; ...
	'yds', @sino_geom_yds, '()'; ...
	'rfov', @sino_geom_rfov, '()'; ...
	};

switch type
case 'fan'
	st = sino_geom_fan(st, varargin{:});
	meth = [meth; { ... % fan-beam methods
		'orbit_short', @sino_geom_orbit_short, '()'; ...
		'gamma', @sino_geom_gamma, '()'; ...
		'gamma_max', @sino_geom_gamma_max, '()'; ...
	}];
case 'par'
	st = sino_geom_par(st, varargin{:});
case 'moj'
	st = sino_geom_moj(st, varargin{:});
otherwise
	fail(['unknown sinotype ' type])
end

if isempty(st.na), st.na = 2 * floor(st.nb * pi/2 / 2); end

st = strum(st, meth);

if st.down ~= 1
	down = st.down; st.down = 1; % trick
	st = st.downsample(down);
end

if streq(type, 'fan') && streq(st.orbit, 'short')
	st.orbit = st.orbit_short;
end

if isempty(st.strip_width) || st.strip_width == 0
	warn 'strip_width=0 is rarely realistic! are you sure?'
end


% sino_geom_orbit_short()
function os = sino_geom_orbit_short(st)
os = 180 + 2*rad2deg(st.gamma_max);


% sino_geom_d()
% sample spacing (radial)
function d = sino_geom_d(st, varargin)
switch st.type
case 'fan'
	if length(varargin), error '?', end
	d = st.ds;
case 'par'
	if length(varargin), error '?', end
	d = st.dr;
case 'moj'
	d = st.dx;
	if length(varargin) % trick: allow st.d(ia) for mojette
		ang = st.ar(varargin{:});
		d = st.dx * max(abs(cos(ang)), abs(sin(ang)));
	end
otherwise
	fail('unknown sinotype "%s"', type)
end


% sino_geom_offset()
% sample offset
function offset = sino_geom_offset(st)
switch st.type
case 'fan'
	offset = st.offset_s;
case 'par'
	offset = st.offset_r;
case 'moj'
	offset = st.offset_r;
otherwise
	fail('unknown sinotype "%s"', type)
end


% sino_geom_dim()
function dim = sino_geom_dim(st)
dim = [st.nb, st.na];
if isempty(st.nb) || isempty(st.na), error 'dim requested without nb,na', end


% sino_geom_ones()
% sinogram of all ones
function out = sino_geom_ones(st)
out = ones(st.dim);


% sino_geom_unitv()
% sinogram with a single ray
function out = sino_geom_unitv(st, ib, ia)
if ~isvar('ib'), ib = round(st.nb/2+1); end
if ~isvar('ia'), ia = round(st.na/2+1); end
out = st.zeros;
out(ib,ia) = 1;


% sino_geom_zeros()
% sinogram of all zeros
function out = sino_geom_zeros(st)
out = zeros(st.dim);


% sino_geom_w()
% 'middle' sample position
function w = sino_geom_w(st)
w = (st.nb-1)/2 + st.offset;


% sino_geom_s()
% sample locations ('radial')
function s = sino_geom_s(st)
s = st.d * ([0:st.nb-1]' - st.w);


% sino_geom_gamma()
% gamma sample values
function gamma = sino_geom_gamma(st, varargin)
switch st.dfs
case 0
	gamma = st.s / st.dsd; % 3rd gen: equiangular
case inf
	gamma = atan(st.s / st.dsd); % flat
otherwise
	error 'not done'
end
gamma = gamma(varargin{:});


% sino_geom_gamma_max()
function gamma_max = sino_geom_gamma_max(st)
gamma_max = max(st.gamma);


% sino_geom_ad()
% angular sample locations (degrees)
function ang = sino_geom_ad(st, varargin)
ang = [0:st.na-1]'/st.na * st.orbit + st.orbit_start;
ang = ang(varargin{:});


% sino_geom_ar()
% angular sample locations (radians)
function ang = sino_geom_ar(st, varargin)
ang = deg2rad(sino_geom_ad(st));
ang = ang(varargin{:});


% sino_geom_downsample()
% down-sample (for testing)
function st = sino_geom_downsample(st, down)
st.down = st.down * down;
st.nb = 2 * round(st.nb / down / 2); % keep it even
st.na = round(st.na / down);

st.strip_width = st.strip_width * down;

switch st.type
case 'fan'
	st.ds = st.ds * down;
case 'par'
	st.dr = st.dr * down;
case 'moj'
	st.dx = st.dx * down;
otherwise
	fail('unknown sinotype "%s"', type)
end


% sino_geom_shape()
% reshape into sinogram array
function sino = sino_geom_shape(st, sino)
sino = reshape(sino, st.nb, st.na, []);


%
% sino_geom_fan()
%
function st = sino_geom_fan(st, varargin);

% defaults
st.orbit = 360; % [degrees]
st.ds		= 1;
st.offset_s	= 0;
st.strip_width	= [];

st.dsd = [];	% dis_src_det
st.dso = [];	% dis_src_iso
st.dod = [];	% dis_iso_det
st.dfs = 0;	% dis_foc_src (3rd gen CT)

subs = { ...
	'ray_spacing', 'ds';
	'channel_offset', 'offset_s';
	'src_det_dis', 'dsd';
	'dis_src_det', 'dsd';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	'ns', 'nb';
	'nbeta', 'na';
	'Dsd', 'dsd';
	'Dso', 'dso';
	'Dod', 'dod';
	'Dfs', 'dfs';
	'obj2det_x', 'dod';
	'obj2det_y', 'dod'
	};
st = vararg_pair(st, varargin, 'subs', subs);

% work out distances
if isempty(st.dsd) + isempty(st.dso) + isempty(st.dod) > 1
	error 'must provide at least two of dsd, dso, dod'
end
if isempty(st.dsd), st.dsd = st.dso + st.dod; end
if isempty(st.dso), st.dso = st.dsd - st.dod; end
if isempty(st.dod), st.dod = st.dsd - st.dso; end
if st.dso + st.dod ~= st.dsd
	error 'bad fan-beam distances'
end

if streq(st.strip_width, 'd') || streq(st.strip_width, 'ds')
	st.strip_width = st.ds;
end


%
% sino_geom_par()
%
function st = sino_geom_par(st, varargin);

% defaults
st.orbit = 180; % [degrees]
st.dr		= 1;
st.offset_r	= 0;
st.strip_width	= [];

subs = { ...
	'nr', 'nb';
	'nphi', 'na';
	'ray_spacing', 'dr';
	'channel_offset', 'offset_r';
	};
st = vararg_pair(st, varargin, 'subs', subs);

if streq(st.strip_width, 'd') || streq(st.strip_width, 'dr')
	st.strip_width = st.dr;
end


% sino_geom_moj()
function st = sino_geom_moj(st, varargin);

% defaults
st.orbit = 180; % [degrees]
st.dx		= 1;
st.offset_r	= 0;
st.strip_width	= []; % ignored

subs = { ...
	'nr', 'nb';
	'nphi', 'na';
	'channel_offset', 'offset_r';
	};
st = vararg_pair(st, varargin, 'subs', subs);


% sino_geom_xds()
% center positions of detectors (for beta = 0)
function xds = sino_geom_xds(st, varargin)
switch st.type
case 'par'
	xds = st.s;
case 'fan'
	switch st.dfs
	case 0 % arc
		gam = st.gamma;
		xds = st.dsd * sin(gam);
	case inf % flat
		xds = st.s;
	otherwise
		error 'not done'
	end
otherwise
	error 'bug'
end
xds = xds + st.source_offset;


% sino_geom_yds()
% center positions of detectors (for beta = 0)
function yds = sino_geom_yds(st, varargin)
switch st.type
case 'par'
	yds = zeros(size(st.s));

case 'fan'
	switch st.dfs
	case 0 % arc
		gam = st.gamma;
		yds = st.dso - st.dsd * cos(gam);
	case inf % flat
		yds = -st.dod * ones(size(st.s));
	otherwise
		error 'not done'
	end
otherwise
	error 'bug'
end


%
% sino_geom_plot()
% picture of the source position / detector geometry
%
function out = sino_geom_plot(st, ig)
switch st.type

case 'par'
	if isvar('ig') && ~isempty(ig)
		im(ig.x, ig.y, ig.mask(:,:,1))
		hold on
	end
	t = linspace(0,2*pi,1001);
	rmax = max(st.s);
	plot(0, 0, '.', rmax * cos(t), rmax * sin(t), '-') % fov circle
	if isvar('ig') && ~isempty(ig)
		hold off
	end

case 'fan'
	x0 = 0;
	y0 = st.dso;
	t = linspace(0,2*pi,100);
	rot = st.ar(1);
	rot = [cos(rot) -sin(rot); sin(rot) cos(rot)];
	p0 = rot * [x0; y0];
	pd = rot * [st.xds'; st.yds'];

	tmp = st.ar + pi/2; % trick: angle beta defined ccw from y axis
	plot(	p0(1), p0(2), 'ys', ... % source
		st.dso * cos(t), st.dso * sin(t), 'c--', ... % source circle
		st.dso * cos(tmp), st.dso * sin(tmp), 'c.', ... % source
		pd(1,:), pd(2,:), 'y.')
	hold on
	if isvar('ig') && ~isempty(ig)
		xmin = min(ig.x); xmax = max(ig.x);
		ymin = min(ig.y); ymax = max(ig.y);
		im(ig.x, ig.y, ig.mask(:,:,1))

		plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	end
	plot(	0, 0, '.', ...
		[pd(1,1) p0(1) pd(1,end)], [pd(2,1) p0(2) pd(2,end)], 'r-', ...
		st.rfov * cos(t), st.rfov * sin(t), 'm:') % fov circle
	xlabel 'x', ylabel 'y', titlef('fov = %g', st.rfov)
	axis equal, axis tight
	hold off

case 'moj'
	if isvar('ig') && ~isempty(ig)
		im(ig.x, ig.y, ig.mask(:,:,1))
		hold on
	end
	t = linspace(0,2*pi,1001);
	rmax = max(st.s);
	rphi = st.nb/2 * st.dx ./ (max(abs(cos(t)), abs(sin(t))));
	plot(0, 0, '.', rmax * cos(t), rmax * sin(t), '-') % fov circle
	plot(0, 0, '.', rphi .* cos(t), rphi .* sin(t), '-m') % fov circle
	if isvar('ig') && ~isempty(ig)
		hold off
	end
	axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)

otherwise
	error 'bug'
end
out = [];


%
% sino_geom_rfov()
% radial FOV
%
function rfov = sino_geom_rfov(st)
switch st.type
case 'par'
	rfov = max(abs(st.s));
case 'fan'
	rfov = st.dso * sin(max(abs(st.gamma)));
otherwise
	fail('not done: %s', st.type)
end
out = [];


%
% sino_geom_taufun()
% projected s/ds, useful for footprint center and support
%
function tau = sino_geom_taufun(st, x, y)
idim = size(x);
jf_equal(idim, size(y))
x = x(:);
y = y(:);
switch st.type
case 'par'
	tau = (x * cos(st.ar) + y * cos(st.ar)) / st.dr;
case 'fan'
	b = st.ar'; % row vector, for outer-product
	xb = x * cos(b) + y * sin(b);
	yb = -x * sin(b) + y * cos(b);
	tangam = (xb - st.source_offset) ./ (st.dso - yb); % e,tomo,fan,L,gam
	switch st.dfs
	case 0 % arc
		tau = st.dsd / st.ds * atan(tangam);
	case inf % flat
		tau = st.dsd / st.ds * tangam;
	otherwise
		fail('dfs %g not done', st.dfs)
	end
otherwise
	fail('not done: %s', st.type)
end


% sino_geom_ge1()
% sinogram geometry for GE lightspeed system
% these numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
%
function geom = sino_geom_ge1(type, varargin)
if ~streq(type, 'ge1'), error 'bug', end

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

units = 'mm'; % default units is mm
scale = 1;
for ii=1:2:length(varargin)
	if streq(varargin{ii}, 'units')
		units = varargin{ii+1};
	end
end
switch units
case 'cm'
	scale = 10;
case {'mm', ''}
otherwise
	fail('unknown units: %s', units)
end
geom = sino_geom('fan', 'nb', 888, 'na', na, ...
	'ds', 1.0239/scale, 'offset_s', 1.25, ...
	'dsd', 949.075/scale, 'dod', 408.075/scale, 'dfs', 0, ...
	'strip_width', 'd', ...
	'units', units, ...
	varargin{:});


% sino_geom_test()
function sino_geom_test
for dfs = [0 inf] % arc flat
	st = sino_geom('fan', 'nb', 888, 'na', 984, ...
		'orbit_start', 20, ...
		'strip_width', 'd', ...
		'source_offset', 0, ...
		'dsd', 949.075, 'dod', 408.075, 'dfs', dfs);
	if im
		clf, st.plot;
	prompt
	end
end
st.ad(2);
st = sino_geom('par', 'strip_width', 'd');
st = sino_geom('moj', 'strip_width', 'd');
st = sino_geom('ge1', 'strip_width', 'd');
st.rfov;
st.downsample(2);
