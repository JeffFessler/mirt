 function img = feldkamp_wt(proj, window, mask, args, varargin)
%function img = feldkamp_wt(proj, window, mask, args, [options])
% modified version by Ajay Paidi
%
% FBP reconstruction of cone-beam tomography data collected with flat detector.
% See feldkamp_example.m for example.
%
% in:
%	proj	[nh,nv,na]	cone-beam projectons (line integrals)
%	window	[npad]		'ramp' (default), or 'hann', or array.
%				if array, then use samples [-K/2, K/2).
%	mask	[nx,ny,nz]	logical array of desired fov support
%	args	char arry	geometry description name/value pairs
% options:
%	'ia_skip' [int]		downsample in angle to save time for tests
% out:
%	img	[nx,ny,nz]		reconstructed image
%
% source_offset has the same units (e.g., mm) as pixel_size etc.
% offset_det_h is an integer or fraction thereof, e.g., 0.25,
% relative to centerline between two central channels.
%
% References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
% Fessler tomography chapter.
%
% Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if nargin == 1 & streq(proj, 'test'), feldkamp_example, return, end
if nargin < 3, help(mfilename), error(mfilename), end
if ~isvar('window') | isempty(window)
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

% extract arguments from name/value pairs
arg_num = @(args,arg) str2num(arg_get(args, arg));
arg_num_def = @(args,arg,def) str2num(arg_get(args, arg, def));

nx = arg_num(args,'nx');
ny = arg_num(args,'ny');
nz = arg_num(args,'nz');
nh = arg_num(args,'nh');
nv = arg_num(args,'nv');
na = arg_num(args,'na');

orbit		= arg_num(args, 'orbit')	* pi / 180;	% to radians
orbit_start	= arg_num(args, 'orbit_start')	* pi / 180;	% to radians
pixel_size	= arg_num(args, 'pixel_size');
ds		= arg_num(args, 'ray_spacing'); % sample spacing [distance]
dt = ds;
strip_width	= arg_num(args, 'strip_width');
center_x	= arg_num_def(args, 'center_x', '0');	% center offset in pixel
center_y	= arg_num_def(args, 'center_y', '0');
center_z	= arg_num_def(args, 'center_z', '0');
offset_source	= arg_num(args, 'offset_source'); % r_off
offset_det_h	= arg_num(args, 'offset_det_h');
offset_det_v	= arg_num(args, 'offset_det_v');
flip_y		= arg_num_def(args, 'flip_y', '1');
Dc		= arg_num(args, 'dis_src_det');
Dd		= arg_num(args, 'dis_iso_det');
Df		= arg_num(args, 'dis_foc_src');
Ds = Dc - Dd; 	% src to "isocenter" distance

if offset_source, error 'only offset_source=0 implemented', end

%
% step 1: weight the projections as in fan-beam case
%
ss = ([-(nh-1)/2:(nh-1)/2]' - offset_det_h) * ds;
tt = ([-(nv-1)/2:(nv-1)/2]' - offset_det_v) * dt;

[ss tt] = ndgrid(ss, tt);
ww1 = Ds * sqrt(1 + (tt/Dc).^2) ./ sqrt(Dc^2 + ss.^2 + tt.^2);
for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia).* ww1;
end
clear ss tt ww1


% step 2: filter the (zero padded) projections
%

npadh = 2^ceil(log2(2*nh-1));
printf('nh=%d npadh=%d', nh, npadh)
projpad = [proj; zeros(npadh-nh,nv,na)]; % padded projections

H = fan_arc_filter(npadh, ds, Dc, window);	% [nb,1]
H = ds * H; % differential for discrete-space convolution vs integral
proj = ifft( fft(projpad) .* repmat(H, [1 nv na]));
proj = proj(1:nh,:,:);
proj = reale(proj);
clear H npad projpad

proj = [proj; zeros(2, nv, na)]; % trick: zero at end saves indexing within loop


%
% step 3: cone-beam backprojection of the filtered views
%

% precompute as much as possible
wx = (nx+1)/2 - center_x;
wy = (ny+1)/2 - center_y;
wz = (nz+1)/2 - center_z;
[xc yc] = ndgrid(([1:nx]-wx) * pixel_size, -flip_y * ([1:ny]-wy) * pixel_size);
zc = ([1:nz]-wz) * pixel_size;
rr = sqrt(xc.^2 + yc.^2);	% [nx,ny]
smax = ((nh-1)/2-abs(offset_det_v)) * ds; % maximum flat panel detector length

gamma_max = atan(smax/Dc);
rmax = Ds * sin(gamma_max);

for myz = 1:nz
	mask(:,:,myz) = mask(:,:,myz) & (rr < rmax);
end
clear wx wy wz rr smax rmax

ia = [0:na-1]';
iz = [0:nz-1]';
betas = orbit_start + orbit * ia / na;	% [na] source angles

% loop over slices
img = zeros(size(mask));
for iz=1:nz
	ticker(mfilename, iz, nz)

	xc2 = xc(mask(:,:,iz));	% [np] pixels within mask for this slice
	yc2 = yc(mask(:,:,iz));

	% loop over each projection angle
	img2 = 0;
	for ia=1:na
		beta = betas(ia);

		x_beta = +xc2 * cos(beta) + yc2 * sin(beta);
		y_beta = -xc2 * sin(beta) + yc2 * cos(beta);

		% detector indices
		mag = Dc ./ (Ds - y_beta);
        weigth = ((Ds + ((abs(iz-(nz/2))/(nz/2))).*y_beta).^4)./((Ds + ((abs(iz-(nz/2))/(nz/2))).*y_beta).^4 + (Ds - ((abs(iz-(nz/2))/(nz/2))).*y_beta).^4);
        %weigth = ((Ds + ((abs(iz-(nz/2))/(nz/2)).^7).*y_beta).^2)./((Ds + ((abs(iz-(nz/2))/(nz/2)).^7).*y_beta).^2 + (Ds - ((abs(iz-(nz/2))/(nz/2)).^7).*y_beta).^2);         
        %weigth = ((Ds + (abs(iz-(nz/2))).*y_beta).^4)./((Ds + ((abs(iz-(nz/2))).*y_beta).^4 + (Ds - ((abs(iz-(nz/2))).*y_beta).^4);
        sprime = mag .* x_beta;
		tprime = mag * zc(iz);
		bh = sprime / ds + (nh+1)/2 + offset_det_h;
		bv = tprime / dt + (nv+1)/2 + offset_det_v;

		% bi-linear interpolation:
		ih = floor(bh); % left bin
		iv = floor(bv);
		igood = 1<=ih & ih<nh & 1<= iv & iv<nv;

		ih(~igood) = nh+1; % trick! point at harmless zeros
		iv(~igood) = 1;

		wr = bh - ih;	% left weight
		wl = 1 - wr;	% right weight
		wu = bv - iv;	% upper weight
		wd = 1 - wu;	% lower weight

		sdim = size(proj);
		ia1 = ia * ones(size(ih));
		p1 =	wl .* proj(sub2ind(sdim, ih,iv,ia1)) + ...
			wr .* proj(sub2ind(sdim, ih+1,iv,ia1));
		p2 =	wl .* proj(sub2ind(sdim, ih,iv+1,ia1)) + ...
			wr .* proj(sub2ind(sdim, ih+1,iv+1,ia1));

		p0 = wu .* p1 + wd .* p2; % vertical interpolation

		p0 = p0 .* mag.^2.* (weigth); % back-projection weighting;
		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = ( orbit / (na/ia_skip)) * embed(img2, mask(:,:,iz));
    
end % iz


%
% apodized filter frequency response
%
function H = fan_arc_filter(n, ds, Dc, window)

h = fan_arc_ramp(n, ds, Dc);
H = reale(fft(fftshift(h)));

if ischar(window)
	if streq(window, 'ramp')
		window = ones(n,1);
	elseif streq(window, 'hann')
		window = hann(n, 'periodic');
	else
		error 'unknown window'
	end
elseif length(window) ~= n
	error 'bad window length'
end

H = H .* fftshift(window);


%
% 'ramp-like' filter for flat fan beam
%
function h = fan_arc_ramp(n, ds, Dc)
n = [-(n/2):(n/2-1)]';
h = zeros(size(n));
h(n==0) = 1 / (4 * ds^2);
odd = mod(n,2) == 1;
h(odd) = -1 ./ (pi * n(odd) * ds).^2;
