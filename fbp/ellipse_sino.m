  function [sino, pos, ang] = ellipse_sino(sg, ells, varargin)
%|function [sino, pos, ang] = ellipse_sino(sg, ells, [options])
%|
%| Create sinogram projection of one or more ellipses.
%| Works for both parallel-beam geometry and for fan-beam geometry.
%|
%| in
%|	sg			sinogram geometry object from sino_geom()
%|	ells	[ne 6]		[centx centy radx rady angle_degrees amplitude]
%|				(or strum motion object; see ellipse_motion.m)
%|
%| options
%|	'oversample'		oversampling factor for emulating "strips"
%|				default 1: just 1 ray per detector element
%|	'xscale'		use -1 to flip in x (not recommended)
%|	'yscale'		use -1 to flip in y (not recommended)
%|
%| out
%|	sino	[nb na]		sinogram
%|	pos	[nb 1]		radial sample positions (or [nb na] if mojette)
%|	ang	[na]		angular sample locations
%|
%| Copyright 2003-10-22, Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if streq(sg, 'test')
	ellipse_sino_test
	prompt
	ellipse_sino_mojette_test
return
end

if isnumeric(sg), ellipse_sino_old(sg, ells, varargin{:}); end % old way

% defaults
arg.xscale = 1;
arg.yscale = 1;
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

switch sg.type
case 'fan'
	[sino pos ang] = ellipse_sino_go(ells, [], [], ...
		sg.nb, sg.ds, sg.offset_s, ...
		sg.na, sg.orbit, sg.orbit_start, ...
		sg.dso, sg.dod, sg.dfs, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, 0);

case 'par'
	[sino pos ang] = ellipse_sino_go(ells, [], [], ...
		sg.nb, sg.dr, sg.offset_r, ...
		sg.na, sg.orbit, sg.orbit_start, ...
		inf, 1, 0, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, 0);

case 'moj'
	[sino pos ang] = ellipse_sino_go(ells, [], [], ...
		sg.nb, [], sg.offset_r, ...
		sg.na, sg.orbit, sg.orbit_start, ...
		inf, 1, 0, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, sg.dx);

otherwise
	fail('sino geom "%s" not done', sg.type)
end


%
% ellipse_sino_go()
%
function [sino, pos, ang] = ellipse_sino_go(ells, ...
	pos, ang, ...
	nb, ds, offset_s, ...
	na, orbit, orbit_start, ...
	dso, dod, dfs, ...
	source_offset, ...
	xscale, yscale, ...
	oversample, mojette);

if isempty(ang)
	ang = deg2rad(orbit_start + [0:na-1]/na * orbit);
end

[pos pos2] = ellipse_sino_pos(pos(:), nb, ds, offset_s, ...
		oversample, mojette, ang);

sino = ellipse_sino_do(ells, pos2, ang(:)', ...
	xscale, yscale, ...
	dso, dod, dfs, source_offset);

if oversample
	sino = downsample2(sino, [oversample 1]);
end


%
% ellipse_sino_pos()
% determine usual and fine "radial" sample positions
%
function [pos_coarse pos_fine] = ...
	ellipse_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)

if isempty(pos)
	if isempty(nb), warning 'nb required when no pos provided', end
	wb = (nb-1)/2 + offset_s;
else
	if ~isempty(nb), warning 'nb ignored when pos provided', end
end

if mojette ~= 0 % tricky mojette radial sampling
	% trick: ray_spacing aka ds comes from dx which is in mojette
	if ~isempty(ds), error 'ds must be empty for mojette case', end
	dt = abs(mojette) * max(abs(cos(ang)), abs(sin(ang)))'; % [1,na]

	if ~isempty(pos), error 'mojette requires empty "pos"', end

	na = length(ang);
	pos_coarse = ([0:nb-1]-wb)' * dt(:)'; % [nb na]
	if nover > 1
		pos_fine = zeros(nover*nb,na);
		for ia=1:na
		tmp = [-(nover-1):2:(nover-1)]' / (2*nover) * dt(ia); % fixed 2008-9-4
			pos_fine(:,ia) = col(outer_sum(tmp, pos_coarse(:,ia)'));
		end
	else
		pos_fine = pos_coarse;
	end

else % ordinary non-mojette sampling

	if isempty(pos)
		pos_coarse = ([0:nb-1]'-wb) * ds; % [nb 1]
	else
		pos_coarse = pos(:); % [nb 1]
		ds = pos(2) - pos(1);
		if any(abs(diff(pos) / ds - 1) > 1e-10)
			error 'uniform spacing required'
		end
	end

	if nover > 1
		% determine fine sampling positions
		% todo: allow different case for trapezoidal rule!
		pos_fine = [-(nover-1):2:(nover-1)]' / (2*nover) * ds; % fixed 2008-9-4
		pos_fine = outer_sum(pos_fine, pos_coarse(:)'); % [nover nb]
		pos_fine = pos_fine(:); % [nb*nover 1]
	else
		pos_fine = pos_coarse;
	end
end


%
% ellipse_sino_do()
% analytical line-integral projections of ellipse (fan-beam in general)
%
function sino = ellipse_sino_do(ells, pos, ang, ...
	xscale, yscale, dso, dod, dfs, source_offset)

nb = size(pos,1);
na = length(ang);

%
% effective radial and angular sample locations in parallel geometry!
%
if isinf(dso)
	if size(pos,2) > 1 % mojette
		rads = pos;
		angs = repmat(ang, [nb 1]); % [nb na]
	else
		[rads angs] = ndgrid(pos, ang); % [nb na]
	end
else % fan
	if size(pos,2) > 1, error 'mojette fan not supported', end
	dis_src_det = dso + dod;

	if isinf(dfs)	% flat detector
		gam = atan(pos / dis_src_det); % "gamma"
	else			% arc detector
		dis_foc_det = dfs + dis_src_det;
		alf = pos / dis_foc_det;
		gam = atan2(	dis_foc_det * sin(alf), ...
				dis_foc_det * cos(alf) - dfs );	% gamma
	end

	rad = dso * sin(gam) + source_offset * cos(gam);	
	rads = rad(:,ones(1,na)); % [nb na]
	angs = outer_sum(gam, ang); % [nb na] gamma + beta
end
clear alf gam rad pos ang

sino = zeros(nb,na);

cangs = cos(angs);
sangs = sin(angs);

% loop over ellipses
ticker reset
if isa(ells, 'strum')
	ne = ells.ne;
else
	ne = size(ells,1);
	if size(ells, 2) ~= 6, fail '6 parameters per ellipse', end
end
for ie = 1:ne
	ticker(mfilename, ie, ne)
	if isa(ells, 'strum')
		ell = ells.ell(ie, na);	% [na 6]
		ell = reshape(ell, [1 na 6]);
		ell = repmat(ell, [nb 1 1]); % [nb na 6]
	else
		ell = ells(ie,:); % [1 6]
		ell = reshape(ell, [1 1 6]);
	end

	cx = xscale * ell(:,:,1);	rx = ell(:,:,3);
	cy = yscale * ell(:,:,2);	ry = ell(:,:,4);
	eang = deg2rad(ell(:,:,5));
	val = ell(:,:,6);

	if yscale == -1
		eang = -eang;
	end
	if xscale == -1
		eang = pi - eang;
	end
	scale = 2 * val .* rx .* ry;

	% square of projected radius:
%	rp2 =	(rx * cos(angs - eang)).^2 + ...
%		(ry * sin(angs - eang)).^2;
	rp2 =	(rx .* (cangs .* cos(eang) + sangs .* sin(eang))).^2 + ...
		(ry .* (sangs .* cos(eang) - cangs .* sin(eang))).^2;
	sp = cx .* cangs + cy .* sangs; % radial shift
	dis2 = (rads - sp).^2; % square of distances from center
	sino = sino + scale ./ rp2 .* sqrt(max(rp2 - dis2,0));
end


 function [sino, pos, ang] = ellipse_sino_old(ells, varargin)
%function [sino, pos, ang] = ellipse_sino_old(ells, sinogeom, [options])
%
% Create sinogram projection of one or more ellipses.
% Works for both parallel-beam geometry and for fan-beam geometry,
% i.e., 3rd generation x-ray ct systems, if the distances are specified.
%
% in:
%	ells	[ne,6]		[centx centy radx rady angle_degrees amplitude]
%	sinogeom struct		sinogram geometry from sino_geom()
%
% options:
%	'oversample'		oversampling factor for emulating "strips"
%	'xscale'		use -1 to flip in x (not recommended)
%	'yscale'		use -1 to flip in y (not recommended)
%
% options for obsolete usage:
%	'pos'	[nb,1]		projection sample locations, typically in mm
%				(arc length along detector)
%				if [], default formed from ds, nb, offset_s,
%					else those three are ignored.
%	'nb'			number of "radial" samples
%	'ds' | 'dr' | 'ray_spacing'	sample spacing, typically in mm
%	'offset_s'		channel offset [unitless]
%
%	'ang'	[na,1]		source angles, ccw from vertical.  aka "beta"
%				in radians.  if [], default formed from
%				na, orbit, orbit_start, else those are ignored
%	'orbit'			[degrees]
%	'orbit_start'
%	'na'
%
%	'dso' | 'Dso' | 'dis_src_iso'	distance from source to isocenter
%					default: Inf, gives parallel-beam
%	'dod' | 'Dod' | 'dis_iso_det'	distance from isocenter to detector
%	'dfs' | 'Dfs' | 'dis_foc_src'	distance from arc focal point to source
%				default is 0, i.e., 3rd generation x-ray CT
%				use Inf for flat fan-beam detector
%	'mojette'	?	if nonzero, use mojette radial sampling.
%				(not recommended except for JF debugging)
%				trick: the nonzero value represents 'dx'
%
% out:
%	sino	[nb,na]		sinogram
%	pos	[nb,1]		radial sample positions (or [nb,na] if mojette)
%	ang	[na]		angular sample locations
%
% Copyright 2003-10-22, Jeff Fessler, The University of Michigan

warn 'obsolete usage; use sino_geom instead; see help ellipse_sino'

% defaults (parallel beam)
arg.pos = [];
arg.nb = [];
arg.ds = 1;
arg.offset_s = 0;

arg.ang = [];
arg.orbit = 360;
arg.orbit_start = 0;
arg.na = [];

arg.mojette = 0.;

arg.dso = inf;
arg.dod = 1;
arg.dfs = 0;

% support very old usage:
% pos, ang, dso, dod, dfs, oversample
if isnumeric(varargin{1})
	[arg varargin] = ellipse_sino_orig(arg, varargin{:});
end

subs = {'dr', 'ds'; 'ray_spacing', 'ds'; 'channel_offset', 'offset_s'; ...
	'Dfs', 'dfs';
	'Dod', 'dod'
	'dis_foc_src', 'dfs';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod'
	};
arg = vararg_pair(arg, varargin, 'subs', subs);

[sino pos ang] = ellipse_sino_go(ells, ...
	arg.pos, arg.ang, ...
	arg.nb, arg.ds, arg.offset_s, ...
	arg.na, arg.orbit, arg.orbit_start, ...
	arg.dso, arg.dod, arg.dfs, ...
	0, ... % source_offset
	arg.xscale, arg.yscale, ...
	arg.oversample, arg.mojette);


%
% ellipse_sino_orig()
%
function [arg, varargin] = ellipse_sino_orig(arg, varargin)
narg = length(varargin);
if narg >= 1, arg.pos = varargin{1}; end
if narg >= 2, arg.ang = varargin{2}; end
if narg >= 3, arg.dso = varargin{3}; end
if narg >= 4, arg.dod = varargin{4}; end
if narg >= 5, arg.dfs = varargin{5}; end
if narg >= 6, arg.oversample = varargin{6}; end
varargin = {varargin{7:end}};

% default is parallel beam
if isempty(arg.dso)
	arg.dso = inf; arg.dod = 1;
end
if isempty(arg.dod)
	error 'dod is required'
end
if isempty(arg.dfs)
	arg.dfs = 0;
end

% default is no over-sampling
if isempty(arg.oversample)
	arg.oversample = 1;
end
if length(arg.ang) == 1
	arg.na = arg.ang;
	arg.ang = [];
end


%
% ellipse_sino_test()
% internal test routine: standard sampling
%
function ellipse_sino_test

down = 4;
ig = image_geom('nx', 512, 'ny', 504, 'dx', 1, 'down', down);
ell = [40 70 50 150 20 10];
[xtrue ell] = ellipse_im(ig, ell, 'oversample', 4);

gf = sino_geom('fan', 'nb', 888, 'na', 984, 'ds', 1.0, 'offset_s', 0.25, ...
	'orbit', 360, 'orbit_start', 0, 'dsd', 949, 'dod', 408, 'down', down);
%	'dfs', inf, 'source_offset', 0.7, ... % flat fan, not working!
gp = sino_geom('par', 'nb', 888, 'na', 984, 'dr', 0.5, 'offset_r', 0.25, ...
	'orbit', gf.orbit, 'orbit_start', gf.orbit_start, 'down', down);

oversample = 8;
sino_mf = ellipse_sino(gf, ell, 'oversample', oversample); % fan
sino_mp = ellipse_sino(gp, ell, 'oversample', 1); % parallel

im pl 3 3
im(1, ig.x, ig.y, xtrue), cbar
%subplot(3), plot(pos, sum(sinof,2))
im(4, gp.s, gp.ad, sino_mp, 'matlab analyt. parallel'), cbar
im(7, gf.s, gf.ad, sino_mf, 'matlab analyt. fan'), cbar

if ~has_aspire, return, end

% compare to aspire analytical parallel beam
t = sprintf('%g %g 0 %g %g %g %g', ell);
f.out = [test_dir 't0'];
com = sprintf('op ellproj,tran %s %d %d 1 %g 0 %g %g %g %s', ...
	f.out, gp.nb, gp.na, gp.dr, gp.offset_r, gp.orbit, 0, t);
os_run(com)
sino_ap = fld_read(f.out); % aspire version of parallel beam

max_percent_diff(sino_ap, sino_mp)

% compare to aspire fan system 13/14
if has_mex_jf

	if oversample > 1, strip_width = gf.ds; else, strip_width = 0; end
	Gf = Gtomo2_dscmex(gf, ig, 'pairs', {'strip_width', strip_width});
	sino_df = Gf * xtrue;
	Gp = Gtomo2_dscmex(gp, ig, 'pairs', {'strip_width', strip_width});
	sino_dp = Gp * xtrue;

	im(5, gp.s, gp.ad, sino_dp, 'aspire sys02 par'), cbar
	im(8, gf.s, gf.ad, sino_df, 'aspire sys14 fan'), cbar
	t = sprintf('par err %g%%', max_percent_diff(sino_mp, sino_dp));
	im(6, gp.s, gp.ad, sino_mp-sino_dp, t), cbar
	t = sprintf('fan err %g%%', max_percent_diff(sino_mf, sino_df));
	im(9, gf.s, gf.ad, sino_mf-sino_df, t), cbar
end


%
% internal test routine: mojette sampling
%
function ellipse_sino_mojette_test

down = 2;
ig = image_geom('nx', 512, 'ny', 496, 'dx', 1, 'down', down);

nb = ceil(max(ig.nx,ig.ny)*sqrt(2)/2)*2;
na = 984/down;
dr = ig.dx;
offset_r = 0.25; % quarter detector
gp = sino_geom('par', 'nb', nb, 'na', na, 'dr', dr, 'offset_r', offset_r);
gm = sino_geom('moj', 'nb', nb, 'na', na, 'dx', ig.dx, 'offset_r', offset_r);

ell = [40 25 150 60 30 10];

oversample = 8;
[sinop pos] = ellipse_sino(gp, ell, 'oversample', oversample);
[sinom1 pos] = ellipse_sino(gm, ell, 'oversample', 1);
[sinom2 pos] = ellipse_sino(gm, ell, 'oversample', oversample);

xtrue = ellipse_im(ig, ell, 'oversample', 4);

im(1, xtrue)
im(4, sinom1, 'mojette 1'), cbar
im(5, sinom2, sprintf('mojette %d', oversample)), cbar
im(6, sinom1-sinom2, 'mojette diff'), cbar
im(7, sinop, 'matlab parallel'), cbar

if has_mex_jf % compare
	if oversample > 1, strip_width = gp.dr; else, strip_width = 0; end
	G = Gtomo2_dscmex(gp, ig, 'pairs', {'strip_width', strip_width});
	sino0 = G * xtrue;

	im(8, sino0, 'aspire sys2'), cbar
	t = sprintf('err %g%%', max_percent_diff(sino0, sinop));
	im(9, sinop-sino0, t), cbar
end
