 function img = feldkamp(proj, window, mask, args, varargin)
%function img = feldkamp(proj, window, mask, args, [options])
%
% FBP reconstruction of cone-beam tomography data collected with
% a circular source trajectory.  
% See feldkamp_example.m for example.
%
% in:
%	proj	[nh,nv,na]	cone-beam projectons (line integrals)
%	window	[npad]		'ramp' (default), or 'hann', or array.
%				if array, then use samples [-K/2, K/2).
%	mask	[nx,ny,nz]	logical array of desired fov support
%	args	char arry	geometry description name/value pairs
%				dis_foc_src: 0 for arc, inf for flat (default)
% options:
%	'ia_skip' [int]		downsample in angle to save time for tests
%
% out:
%	img	[nx,ny,nz]		reconstructed image
%
% offset_source has the same units (e.g., mm) as pixel_size etc.
% offset_det_h is an integer or fraction thereof, e.g., 0.25,
% relative to centerline between two central channels.
%
% References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
% Fessler tomography chapter.
%
% Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, The University of Michigan
% arc detector case contributed by Yingying Zhang 2005-6-13

if nargin == 1 && streq(proj, 'test'), run_mfile_local('feldkamp_example'), return, end
if nargin < 3, help(mfilename), error(mfilename), end
if ~isvar('window') || isempty(window)
	window = 'ramp';
end

arg.ia_skip = 1;
arg = vararg_pair(arg, varargin);

% extract arguments from name/value pairs
arg_num = @(args,arg) str2num(arg_get(args, arg));
arg_num_def = @(args,arg,def) str2num(arg_get(args, arg, def));

nx = arg_num(args, 'nx');
ny = arg_num(args, 'ny');
nz = arg_num_def(args, 'nz', size(proj,2));
nh = arg_num_def(args, 'nh', size(proj,1));
nv = arg_num_def(args, 'nv', size(proj,2));
na = arg_num_def(args, 'na', size(proj,3));

orbit		= arg_num(args, 'orbit')	* pi / 180;	% to radians
orbit_start	= arg_num(args, 'orbit_start')	* pi / 180;	% to radians
pixel_size	= arg_num(args, 'pixel_size');
ds		= arg_num(args, 'ray_spacing'); % sample spacing [distance]
dt = ds;
strip_width	= arg_num(args, 'strip_width');
center_x	= arg_num_def(args, 'center_x', '0');	% center offset in pixel
center_y	= arg_num_def(args, 'center_y', '0');
center_z	= arg_num_def(args, 'center_z', '0');
offset_source	= arg_num_def(args, 'offset_source', 0); % r_off
t		= arg_num_def(args, 'source_offset', 0);
offset_det_h	= arg_num_def(args, 'offset_det_h', t);
offset_det_v	= arg_num_def(args, 'offset_det_v', 0);
flip_y		= arg_num_def(args, 'flip_y', '1');
Dsd		= arg_num(args, 'dis_src_det');
Dod		= arg_num(args, 'dis_iso_det');
Dfs		= arg_num_def(args, 'dis_foc_src', inf); % default flat
Dso = Dsd - Dod; 	% src to "isocenter" distance

if offset_source, error 'only offset_source=0 implemented', end

%
% step 1: weight the projections as in fan-beam case
%
ss = ([-(nh-1)/2:(nh-1)/2]' - offset_det_h) * ds;
tt = ([-(nv-1)/2:(nv-1)/2]' - offset_det_v) * dt;

[ss tt] = ndgrid(ss, tt);
if isinf(Dfs) % flat
	ww1 = Dso * sqrt(1 + (tt/Dsd).^2) ./ sqrt(Dsd^2 + ss.^2 + tt.^2);
elseif Dfs == 0 % arc
	ww1 = (Dso/Dsd) * cos(ss ./ (Dsd * sqrt(1 + (tt/Dsd).^2)));
else
	error('other configurations not implemented.')
end   

for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia) .* ww1;
end
clear ss tt ww1


%
% step 2: filter the (zero padded) projections
%

npadh = 2^ceil(log2(2*nh-1));
printf('nh=%d npadh=%d', nh, npadh)

if isinf(Dfs)
	H = fan_filter('flat', npadh, ds, [], window);	% [nb,1]
elseif Dfs == 0
	H = fan_filter('arc', npadh, ds, Dsd, window);
end
H = ds * H; % differential for discrete-space convolution vs integral

%projpad = [proj; zeros(npadh-nh,nv,na)]; % padded projections
proj = ifft_sym( fft(proj, npadh, 1) .* repmat(H, [1 nv na]) );
%proj = reale(proj);
proj = proj(1:nh,:,:);
clear H npad %projpad

proj = [proj; zeros(2, nv, na)]; % trick: zero at end saves indexing within loop


%
% step 3: cone-beam backprojection of the filtered views
%

% precompute as much as possible
wx = (nx+1)/2 + center_x;
wy = (ny+1)/2 + center_y;
wz = (nz+1)/2 + center_z;
[xc yc] = ndgrid(([1:nx]-wx) * pixel_size, -flip_y * ([1:ny]-wy) * pixel_size);
zc = ([1:nz]-wz) * pixel_size; % right here assumes cubic voxels; could change.
rr = sqrt(xc.^2 + yc.^2);	% [nx,ny]
smax = ((nh-1)/2-abs(offset_det_h)) * ds; % maximum detector s coordinate

if isinf(Dfs)
	gamma_max = atan(smax/Dsd);
elseif Dfs == 0
	gamma_max = smax / Dsd;
end

rmax = Dso * sin(gamma_max);

for myz = 1:nz
	mask(:,:,myz) = mask(:,:,myz) & (rr < rmax);
end
[rmax sum(col(mask(:,:,1)))]
clear wx wy wz rr smax rmax

ia = [0:na-1]';
iz = [0:nz-1]';
betas = orbit_start + orbit * ia / na;	% [na] source angles

% loop over slices
img = zeros(size(mask));
ticker reset
for iz=1:nz

	xc2 = xc(mask(:,:,iz));	% [np] pixels within mask for this slice
	yc2 = yc(mask(:,:,iz));

	% loop over each projection angle
	img2 = 0;
	for ia=1:na
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc2 * cos(beta) + yc2 * sin(beta);
		y_beta = -xc2 * sin(beta) + yc2 * cos(beta);

		% detector indices
		mag = Dsd ./ (Dso - y_beta);
        
		if isinf(Dfs)
			sprime = mag .* x_beta;
		elseif (Dfs == 0)
			d_loop = Dso - y_beta;
			r_loop = x_beta - offset_source;
			sprime = Dsd * atan2(r_loop, d_loop);
		end
        
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
        
		if isinf(Dfs)
			p0 = p0 .* mag.^2; % back-projection weighting for flat
		elseif Dfs == 0
			p0 = p0 .* (Dsd.^2) ./ (r_loop.^2 + d_loop.^2);
		end
        
		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask(:,:,iz));
end % iz

img = (0.5 * orbit / (na/arg.ia_skip)) * img; % final "\der angle" scale


%
% apodized filter frequency response
%
function H = fan_filter(type, n, ds, Dsd, window)

if streq(type, 'flat')
	h = fbp_ramp('flat', n, ds);
else
	h = fbp_ramp('arc', n, ds, Dsd);
end
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
