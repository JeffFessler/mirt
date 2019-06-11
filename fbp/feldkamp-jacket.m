  function img = feldkamp(cg, ig, proj, varargin)
%|function img = feldkamp(cg, ig, proj, varargin)
%|
%| version modified by james balter and others to try accelereys "jacket"
%| UNDER DEVELOPMENT!
%|
%| FBP reconstruction of cone-beam tomography data collected with
%| a circular source trajectory.
%| See feldkamp_example.m for example.
%|
%| in:
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns,nt,na]	cone-beam projection views (line integrals)
%|
%| options:
%|	'window' [npad]		'ramp' (default), or 'hann', or array.
%|				if array, then use samples [-K/2, K/2).
%|	'offset_source'		distance from isocenter to perpendicular ray
%|				[the same units (e.g., mm) as pixel_size etc.]
%|				caution: probably should not be used
%|	'ia_skip' [int]		downsample in angle to save time for tests
%|	'use_mex' 0|1		backprojector: 0 for matlab, 1 for mex (default)
%|	'nthread' 		default: jf('ncore')
%|
%| out:
%|	img	[nx,ny,nz]	reconstructed image
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

if isnumeric(cg)
	img = feldkamp_old(cg, ig, proj, varargin{:});
return
end

% defaults
arg.use_mex = 1;
arg.window = 'ramp';
arg.offset_source = 0; % r_off, distance between rotation iso-center
			% and ray from source that is orthogonal to detector.
arg.ia_skip = 1;
arg.nthread = jf('ncore');
arg = vararg_pair(arg, varargin);

if cg.pitch ~= 0 || any(cg.zshifts ~= 0)
	fail('sorry, helical CT unsupported')
end

img = feldkamp_do(proj, ...
	cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
	cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
	ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
	[ig.offset_x ig.offset_y ig.offset_z], ...
	arg.window, arg.ia_skip, arg.use_mex, arg.nthread);
end % feldkamp()


 function img = feldkamp_old(proj, mask, varargin)
%function img = feldkamp_old(proj, mask, [options])
%
% FBP reconstruction of cone-beam tomography data collected with
% a circular source trajectory.
% See feldkamp_example.m for example.
%
% in:
%	proj	[ns,nt,na]	cone-beam projection views (line integrals)
%	mask	[nx,ny]		logical array of desired (2d) fov support
%
% options:
%	window	[npad]		'ramp' (default), or 'hann', or array.
%				if array, then use samples [-K/2, K/2).
%	'dx' 'dy' 'dz'		voxel dimensions
%	'nz'			# of slices
%	'center_xyz'		image center [pixels] (default: 0)
%	'ds' 'dt'		detector sample spacing (horiz, vert)
%	'offset_st'		detector offset_s,t [integer number of samples
%				or fraction thereof, for 1/4-detector offset.
%	'orbit'			default 360 degrees
%	'orbit_start'		default 0 degrees
%	'dis_src_det'
%	'dis_iso_det'
%	'dis_foc_src'		default inf (flat), or 0 for arc
%	'offset_source'		distance from isocenter to perpendicular ray
%				[the same units (e.g., mm) as pixel_size etc.]
%	'ia_skip' [int]		downsample in angle to save time for tests
%	'use_mex' 0|1		backprojector: 0 for matlab, 1 for mex (default)
%
% out:
%	img	[nx,ny,nz]		reconstructed image
%
% offset_s is an integer or fraction thereof, e.g., 0.25,
% relative to centerline between two central channels.
%
% References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
% Fessler tomography chapter.
%
% Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, University of Michigan
% arc detector case contributed by Yingying Zhang 2005-6-13

warn 'this is obsolete use of feldkamp.m - use new style!'

% defaults
arg.use_mex = 1; % force MATLAB
arg.window = 'ramp';
arg.dx = nan;
arg.dy = [];
arg.dz = [];
arg.nz = size(proj,2);
arg.center_xyz = [0 0 0]; % center_x,y,z (pixels)
arg.flip_y = 1;

arg.dsd = []; % 'dis_src_det'
arg.dod	= []; %	'dis_iso_det'
arg.dfs = inf; % 'dis_foc_src' default flat
arg.offset_source = 0; % r_off, distance between rotation iso-center
			% and ray from source that is orthogonal to detector.

arg.ds = nan;
arg.dt = [];
arg.offset_st = [0 0]; % offset_s,t
arg.ia_skip = 1;
arg.orbit = 360;
arg.orbit_start = 0;
arg = vararg_pair(arg, varargin, 'subs', ...
	{'dis_src_det', 'dsd';
	'dis_iso_det', 'dod'; 
	'dis_foc_src', 'dfs'});
if isnan(arg.ds) || isnan(arg.dx), fail('ds and dx required'), end

if isempty(arg.dy), arg.dy = -arg.dx; end
if isempty(arg.dz), arg.dz = arg.dx; end
if isempty(arg.dt), arg.dt = arg.ds; end

arg.dso = arg.dsd - arg.dod; 	% src to "isocenter" distance

if ~isempty(arg.zshifts), error 'todo: helical feldkamp not done', end

img = feldkamp_do(proj, ...
	arg.ds, arg.dt, arg.offset_st(1), arg.offset_st(2), ...
	arg.offset_source, ...
	arg.dsd, arg.dso, arg.dfs, ...
	arg.orbit, arg.orbit_start, ...
	mask, arg.nz, arg.dx, arg.dy, arg.dz, arg.center_xyz, ...
	arg.window, arg.ia_skip, arg.use_mex, 1);

end % feldkamp_old()


%
% feldkamp_do()
%
function img = feldkamp_do(proj, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	window, ia_skip, use_mex, nthread)

% step 1: weight the projections as in fan-beam case
proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs);

% step 2: filter the (zero padded) projections
[ns nt na] = size(proj);
if use_mex
	pad1 = 0;
else
	pad1 = 2;
end
proj = feldkamp_filter(proj, window, dsd, dfs, ds, pad1);

if pad1 % trick: zero at end saves indexing in loop
%	proj = [proj; zeros(2, nt, na)];
	proj(:,end+1,:) = 0;
end

% step 3: cone-beam backprojection of the filtered views
cpu etic
img = feldkamp_back(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	ia_skip, use_mex, nthread);
cpu etoc 'fdk backprojection cpu time:'

end % feldkamp_do()


%
% feldkamp_weight1()
% step 1: weight the projections as in fan-beam case
%
function proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, ...
	dsd, dso, dfs);
[ns nt na] = size(proj);
ss = ([-(ns-1)/2:(ns-1)/2]' - offset_s) * ds;
tt = ([-(nt-1)/2:(nt-1)/2]' - offset_t) * dt;

[ss tt] = ndgrid(ss, tt);
if isinf(dfs) % flat
	ww1 = dso * sqrt(1 + (tt/dsd).^2) ./ sqrt(dsd^2 + ss.^2 + tt.^2);
elseif dfs == 0 % arc
	ww1 = (dso/dsd) * cos(ss ./ (dsd * sqrt(1 + (tt/dsd).^2)));
else
	error 'other configurations not implemented'
end

for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia) .* ww1;
end
end % feldkamp_weight1()


%
% feldkamp_filter()
% step 2: filter the (zero padded) projections
%
function proj = feldkamp_filter(proj, window, dsd, dfs, ds, pad1)
[ns nt na] = size(proj);
npadh = 2^ceil(log2(2*ns-1));
printf('ns=%d npadh=%d', ns, npadh)

if isinf(dfs)
	H = fan_filter('flat', npadh, ds, [], window);	% [nb,1]
elseif dfs == 0
	H = fan_filter('arc', npadh, ds, dsd, window);
end
H = ds * H; % differential for discrete-space convolution vs integral

%% START JACKET
GPU = 1;
if GPU,
    proj = gsingle( proj );
    H = gsingle( H );
    proj_large = gzeros( npadh, size(proj,2), size(proj,3) );
else,
    proj_large = zeros( npadh, size(proj,2), size(proj,3) );
end
%% END JACKET

tic;
proj_large(1:size(proj,1),:,:) = proj;
proj = ifft( fft(proj_large) .* repmat(H, [1 nt na]) );
proj = proj(1:(ns+pad1),:,:); % trick: extra zero at end saves indexing in loop
toc

%% START JACKET
if GPU,
    proj = double( proj );
end
%% END JACKET

end % feldkamp_filter()


%
% fan_filter()
% apodized filter frequency response
%
function H = fan_filter(type, n, ds, dsd, window)

if streq(type, 'flat')
	h = fbp_ramp('flat', n, ds);
else
	h = fbp_ramp('arc', n, ds, dsd);
end
H = reale(fft(fftshift(h)));

if ischar(window)
	if streq(window, 'ramp')
		window = ones(n,1);
	elseif streq(window, 'hann')
		window = ir_hann_periodic(n);
	else
		error 'unknown window'
	end
elseif length(window) ~= n
	error 'bad window length'
end

H = H .* fftshift(window);
end % fan_filter()


%
% feldkamp_back()
% step 3: cone-beam backprojection of the filtered views
%
function img = feldkamp_back(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	ia_skip, use_mex, nthread)
[nx ny] = size(mask);

betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles

% mex backprojector
if use_mex
	proj = single(proj);
	nthread = int32(nthread);

	if use_mex == 1 % loop in mex
		proj = permute(proj, [2 1 3]); % ts
		tmp = jf_mex('fdk,ts,back', int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
                		dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj, betas, nthread);
		img = double6(tmp);
		img = permute(img, [2 3 1]); % zxy -> xyz

	elseif use_mex == 2 % loop in matlab (for testing)
		proj = permute(proj, [2 1 3]); % ts
		img = 0;
		for ia=1:na
			ticker(mfilename, ia, na)
			% note: 2006-5-30: replaced -dy with dy
			tmp = jf_mex('fdk,ts,back', ...
				int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
                		dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj(:,:,ia), betas(ia), nthread);
			tmp = double6(tmp);
			img = img + tmp;
		end
		img = permute(img, [2 3 1]); % zxy -> xyz

	elseif use_mex == 3 % fdk,st (for testing only - slower!)
		tmp = jf_mex('fdk,st,back', int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
                		dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj, betas, nthread);
		img = double6(tmp);

	else
		fail 'bug'
	end

	% final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
return
end

% matlab backprojector (slower)

fprintf('We are using slower version!\n');

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

% loop over slices
GPU = 0;
if GPU,
    img = gzeros([size(mask) nz]) + i;
    img = img - i;
else,
    img = zeros([size(mask) nz]);
end;
ticker reset
tic;
for iz=1:nz

	% loop over each projection angle
	img2 = 0;
	for ia=1:na
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc * cos(beta) + yc * sin(beta);
		y_beta = -xc * sin(beta) + yc * cos(beta);

		% detector indices
		mag = dsd ./ (dso - y_beta);

		if isinf(dfs)
			sprime = mag .* x_beta;
		elseif (dfs == 0)
			d_loop = dso - y_beta;
			r_loop = x_beta - offset_source;
			sprime = dsd * atan2(r_loop, d_loop);
		end

		tprime = mag * zc(iz);
		bh = sprime / ds + (ns+1)/2 + offset_s;
		bv = tprime / dt + (nt+1)/2 + offset_t;

		% bi-linear interpolation:
		ih = floor(bh); % left bin
		iv = floor(bv);
		igood = 1<=ih & ih<ns & 1<= iv & iv<nt;

		ih(~igood) = ns+1; % trick! point at harmless zeros
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

		if isinf(dfs)
			p0 = p0 .* mag.^2; % back-projection weighting for flat
		elseif dfs == 0
			p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + d_loop.^2);
		end

		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask);
end % iz
toc

img = (0.5 * deg2rad(orbit) / (na/ia_skip)) * img; % final "\der angle" scale
img = double( img );

end % feldkamp_back()
