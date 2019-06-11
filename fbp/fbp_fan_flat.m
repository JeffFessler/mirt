 function img = fbp_fan_flat(sino, G, window, varargin)
%function img = fbp_fan_flat(sino, G, window, varargin)
%|
%| FBP reconstruction of fan-beam tomography data for a flat detector.
%| (Equidistant case.)  Extracts parameters from Gtomo2_dsc object.
%| See fbp_fan_flat_example.m for example.
%|
%| in:
%|	sino	[nb,na]		sinogram (line integrals)
%|	G	[nd,np]		Gtomo2_dscmex object
%|	window	[npad]		'ramp', or 'hann', or array (default: 'ramp')
%|				if array, then use samples [-K/2, K/2)
%|	varargin		options
%|		'ia_skip' [int]	downsample in angle to save time for quick tests
%| out:
%|	img	[nx,ny]		reconstructed image
%|
%| source_offset has the same units (e.g., mm) as pixel_size etc.
%| channel_offset is an integer or fraction thereof, e.g., 0.25.
%|	gamma = dgamma * (ib - (nb-1)/2 - channel_offset)
%| channel offset is relative to centerline between two central channels.
%|
%| References:
%| Gullberg et al, T-MI, Mar. 1986 (gullberg:86:raf).
%| Kak & Slaney, p. 77.  Fessler tomography chapter.
%|
%| Copyright 2004-5 by Patty Laskowsky, Nicole Caparanis, Taka Masuda
%| and Jeff Fessler, The University of Michigan

if nargin == 1 && streq(sino, 'test'), run_mfile_local('fbp_fan_flat_example'), return, end
if nargin < 2, help(mfilename), error(mfilename), end
if ~isvar('window') || isempty(window)
	window = 'ramp';
end

ia_skip = 1;
while(length(varargin))
	arg = varargin{1};
	if streq(arg, 'ia_skip')
		ia_skip = varargin{2};
		varargin = {varargin{3:end}};
	else
		error 'unknown optional argument'
	end
end

% extract arguments from system 14 Gtomo2_dsc object
%nx = G.nx;
%ny = G.ny;
[nx ny] = size(G.arg.mask);
nb = G.arg.nb;
na = G.arg.na;

arg_num = @(G,arg) str2num(arg_get(G.arg.args, arg));
arg_num_def = @(G,arg,def) str2num(arg_get(G.arg.args, arg, def));
orbit		= arg_num(G, 'orbit')		* pi / 180; % to radians
orbit_start	= arg_num(G, 'orbit_start')	* pi / 180; % to radians
pixel_size	= arg_num(G, 'pixel_size');
ds		= arg_num(G, 'ray_spacing'); % sample spacing [distance]
%strip_width	= arg_num(G, 'strip_width');
Dsd		= arg_num(G, 'src_det_dis');	% dis_src_det
Dod		= arg_num(G, 'obj2det_x');	% dis_iso_det
obj2det_y	= arg_num(G, 'obj2det_y');
Dfs		= arg_num_def(G, 'dis_foc_src', Inf);	% Inf for flat det.
center_x	= arg_num_def(G, 'center_x', '0');	% center offset in pixel
center_y	= arg_num_def(G, 'center_y', '0');
source_offset	= arg_num_def(G, 'source_offset', 0); % r_off
channel_offset	= arg_num_def(G, 'channel_offset', 0);
flip_y		= arg_num_def(G, 'flip_y', '1');

if obj2det_y ~= Dod, error 'only circular orbit implemented', end

%
% fan beam FBP step 1: weight the sinogram
%

nn = [-(nb-1)/2:(nb-1)/2]' - channel_offset;
Dso = Dsd - Dod;		% src to "isocenter" distance
ss = ds * nn;
gam = atan(ss / Dsd);
w1 = abs(Dso * cos(gam) - source_offset * sin(gam)) / Dsd;  % 1D weighting
sino = sino .* w1(:,ones(1,na));
clear nn gam w1 ss


%
% fan beam FBP step 2: filter the (zero padded) sinogram
% trick: extra zero column saves linear interpolation indexing within loop!
%
sino = fbp2_sino_filter('flat', sino, 'ds', ds, 'extra', 1, 'window', window);

%
% fan beam FBP step 3: backproject the filtered sinogram
%

% precompute as much as possible
wx = (nx+1)/2 - center_x;
wy = (ny+1)/2 - center_y;
[xc yc] = ndgrid(([1:nx]-wx) * pixel_size, -flip_y * ([1:ny]-wy) * pixel_size);
rr = sqrt(xc.^2 + yc.^2);	% [nx,ny]

smax = ((nb-1)/2-abs(channel_offset)) * ds; %maximum flat panel detector length
gamma_max = atan(smax / Dsd);
rmax = Dso * sin(gamma_max);
mask = G.mask & (rr < rmax);
xc = xc(mask(:));	% [np] pixels within mask
yc = yc(mask(:));
clear wx wy rr smax rmax

ia = [0:na-1]';
betas = orbit_start + orbit * ia / na;	% [na]

img = 0;

% loop over each projection angle
for ia=1:ia_skip:na
	ticker(mfilename, ia, na)

	beta = betas(ia);
	d_loop = Dso + xc * sin(beta) - yc * cos(beta); % Dso - y_beta
	x_beta = xc * cos(beta) + yc * sin(beta);
	w2 = Dsd^2 ./ d_loop.^2;	% [np] back projection weighting  
	sprime_ds = (Dsd / ds) * (x_beta - source_offset) ./ d_loop; % s' / ds
	bb = sprime_ds + ((nb+1)/2 + channel_offset); % [np] bin "index"
 
	% nearest neighbor interpolation:
%	ib = round(bb);
%	if any(ib < 1 | ib > nb), error 'bug', end
%	% trick: make out-of-sinogram indices point to those extra zeros
%	ib(ib < 1 | ib > nb) = nb+1;
%	img = img + sino(ib, ia) ./ L2;

	% linear interpolation:
	il = floor(bb);	% left bin
%	if any(il < 1 | il >= nb), error 'bug', end
	wr = bb - il;	% left weight
	wl = 1 - wr;	% right weight
	img = img + (wl .* sino(il, ia) + wr .* sino(il+1, ia)) .* w2;
end

img = (0.5 * orbit / (na/ia_skip)) * embed(img, mask);
