 function img = fbp_fan_arc(sino, G, window, varargin)
%function img = fbp_fan_arc(sino, G, window, [options])
%|
%| FBP reconstruction of fan-beam tomography data collected with arc detector.
%| Assumes system 14 geometry and extracts parameters from Gtomo2_dsc object.
%| Assumes that the arc is focused on the source, e.g., 3rd generation X-ray CT.
%| See fbp_fan_arc_example.m for example.
%|
%| in:
%|	sino	[nb,na,nz]	sinogram(s) (line integrals)
%|	G	[nd,np]		Gtomo2_dscmex object, or empty if varargin used!
%|	window	[npad]		'ramp', or 'hann', or array (default: 'ramp')
%|				if array, then use samples [-K/2, K/2)
%| options:
%|	'ia_skip' [int]		downsample in angle to save time for quick tests
%|	'rmax' [value]		maximum radius within FOV
%|				(default: usual fully-sampled FOV)
%|	'mask'
%|	'nx'
%|	'ny'
%|	'orbit'
%|	'orbit_start'
%|	'pixel_size'
%|	'ds'
%|	'strip_width'
%|	'center_x'
%|	'center_y'
%|	'source_offset'
%|	'channel_offset'
%|	'flip_x'
%|	'flip_y'
%|
%|	'dsd		= [];	% dis_src_det
%|	'dod		= [];	% dis_iso_det
%|	'obj2det_y	= [];
%|
%| out:
%|	img	[nx,ny,nz]	reconstructed image(s)
%|
%| source_offset has the same units (e.g., mm) as pixel_size etc.
%| channel_offset is an integer or fraction thereof, e.g., 0.25.
%| channel offset is relative to centerline between two central channels.
%|
%| References:
%| Gullberg et al, T-MI, Mar. 1986 (gullberg:86:raf).
%| Kak & Slaney, p. 77. Fessler tomography chapter.
%|
%| Copyright December 2001. Idris Elbakri, The University of Michigan
%| Accelerated by Jeff Fessler 2002-2-19, rewritten 2004-5-21.

if nargin == 1 && streq(sino, 'test'), run_mfile_local('fbp_fan_arc_example'), return, end
if nargin < 2, help(mfilename), error(mfilename), end
if ~isvar('window') || isempty(window)
	arg.window = 'ramp';
else
	arg.window = window;
end

% defaults
arg.mask = [];
arg.nx = 2^9;
arg.ny = 2^9;
arg.ia_skip = 1;
arg.rmax = [];
arg.orbit	= pi;
arg.orbit_start	= 0;
arg.pixel_size	= 1;
arg.ds		= 1;
arg.strip_width	= []; % not needed really
arg.center_x	= 0;
arg.center_y	= 0;
arg.source_offset	= 0;
arg.channel_offset	= 0;
arg.flip_x		= 1;
arg.flip_y		= 1;

arg.dsd		= [];	% dis_src_det
arg.dod		= [];	% dis_iso_det
arg.obj2det_y	= [];

arg = vararg_pair(arg, varargin, 'subs', ...
	{'src_det_dis', 'dsd'; 'Dsd', 'dsd';
	'obj2det_x', 'dod'; 'Dod', 'dod'});

% extract arguments from system 14 Gtomo2_dsc object, if provided
if ~isempty(G)
	arg = fbp_fan_arc_parse_struct(arg, G);
	if size(sino,1) ~= G.arg.nb || size(sino,2) ~= G.arg.na
		error 'bad sino size'
	end
end

if arg.obj2det_y ~= arg.dod, error 'only circular orbit implemented', end

img = fbp_fan_arc_do(sino, arg.window, ...
	arg.mask, arg.nx, arg.ny, ...
	arg.ia_skip, arg.rmax, ...
	arg.orbit, arg.orbit_start, ...
	arg.pixel_size, arg.ds, arg.strip_width, ...
	arg.center_x, arg.center_y, ...
	arg.source_offset, arg.channel_offset, ...
	arg.flip_x, arg.flip_y, ...
	arg.dsd, arg.dod);

function img = fbp_fan_arc_do(sino, window, ...
	mask, nx, ny, ...
	ia_skip, rmax, ...
	orbit, orbit_start, ...
	pixel_size, ds, strip_width, ...
	center_x, center_y, ...
	source_offset, channel_offset, ...
	flip_x, flip_y, ...
	dsd, dod);

[nb na nz] = size(sino);

%
% fan beam FBP step 1: weight the sinogram(s)
%

nn = [-(nb-1)/2:(nb-1)/2]' - channel_offset;
dso = dsd - dod;	% src to "isocenter" distance
ss = nn * ds;		% sample locations
w1 = abs(dso * cos(ss/dsd) - source_offset * sin(ss/dsd)) / dsd; % 1D weighting
sino = sino .* repmat(w1, [1 na nz]);
clear nn w1 ss


%
% fan beam FBP step 2: filter the sinogram
% trick: extra zero column saves linear interpolation indexing within loop!
%

sino = fbp2_sino_filter('arc', sino, 'ds', ds, 'dsd', dsd, ...
	'extra', 1, 'window', window);

%
% fan beam FBP step 3: backproject the filtered sinogram
%

% precompute as much as possible
wx = (nx+1)/2 - center_x;
wy = (ny+1)/2 - center_y;
[xc yc] = ndgrid(flip_x * ([1:nx]-wx) * pixel_size, ...
		-flip_y * ([1:ny]-wy) * pixel_size);
rr = sqrt(xc.^2 + yc.^2);	% [nx,ny]

smax = ((nb-1)/2-abs(channel_offset)) * ds;
if isempty(rmax)
	rmax = dso * sin(smax / dsd); % alpha_max
end
mask = mask & (rr < rmax);
xc = xc(mask(:));	% [np] pixels within mask
yc = yc(mask(:));
clear wx wy rr smax

ia = [0:na-1]';
betas = orbit_start + orbit * ia / na;	% [na]
wb = (nb+1)/2 + channel_offset;

img = 0;

% loop over each projection angle
for ia=1:ia_skip:na
	ticker(mfilename, ia, na)

	beta = betas(ia);
	d_loop = dso + xc * sin(beta) - yc * cos(beta);	% dso - y_beta
	r_loop = xc * cos(beta) + yc * sin(beta) - source_offset; % x_beta-roff
	w2 = dsd^2 ./ (d_loop.^2 + r_loop.^2);	% [np] image-domain weighting
	sprime_ds = (dsd/ds) * atan2(r_loop, d_loop); % s' / ds

	bb = sprime_ds + wb; % [np] bin "index"

	% nearest neighbor interpolation:
%	ib = round(bb);
%	if any(ib < 1 | ib > nb), error 'bug', end
%	% trick: make out-of-sinogram indices point to those extra zeros
%%	ib(ib < 1 | ib > nb) = nb+1;
%	img = img + sino(ib, ia) ./ L2;

	% linear interpolation:
	il = floor(bb);	% left bin
	ir = 1+il;	% right bin

	% deal with truncated sinograms
	ig = il >= 1 & ir <= nb;
	il(~ig) = nb+1;
	ir(~ig) = nb+1;
%	if any(il < 1 | il >= nb), error 'bug', end

	wr = bb - il;	% left weight
	wl = 1 - wr;	% right weight
	if nz > 1
		wr = repmat(wr, [1 nz]);
		wl = repmat(wl, [1 nz]);
		img = img + (wl .* squeeze(sino(il, ia, :)) ...
			+ wr .* squeeze(sino(ir, ia, :))) .* repmat(w2, [1 nz]);
	else
		img = img + (wl .* sino(il, ia) + wr .* sino(ir, ia)) .* w2;
	end
end

img = pi / (na/ia_skip) * embed(img, mask);


% fbp_fan_arc_parse_struct()
function arg = fbp_fan_arc_parse_struct(arg, G);

arg_num = @(G, arg) str2num(arg_get(G.arg.args, arg));
arg_num_def = @(G, arg, def) str2num(arg_get(G.arg.args, arg, def));

arg.orbit	= arg_num(G, 'orbit')		* pi / 180; % to radians
arg.orbit_start	= arg_num(G, 'orbit_start')	* pi / 180; % to radians
arg.pixel_size	= arg_num(G, 'pixel_size');
arg.ds		= arg_num(G, 'ray_spacing'); % sample spacing [distance]
arg.strip_width	= arg_num(G, 'strip_width');
arg.dsd		= arg_num(G, 'src_det_dis');	% dis_src_det
arg.dod		= arg_num(G, 'obj2det_x');	% dis_iso_det
arg.obj2det_y	= arg_num(G, 'obj2det_y');
arg.center_x	= arg_num_def(G, 'center_x', '0'); % center offset in pixel
arg.center_y	= arg_num_def(G, 'center_y', '0');
arg.source_offset	= arg_num(G, 'source_offset'); % r_off
arg.channel_offset	= arg_num(G, 'channel_offset');
arg.flip_x		= arg_num_def(G, 'flip_x', '1');
arg.flip_y		= arg_num_def(G, 'flip_y', '1');

%arg.nx = G.nx;
%arg.ny = G.ny;

arg.mask = G.arg.mask;
[arg.nx arg.ny] = size(arg.mask);
