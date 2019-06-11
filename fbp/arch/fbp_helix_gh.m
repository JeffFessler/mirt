  function [img sino] = fbp_helix_gh(cg, ig, proj, varargin)
%|function [img sino] = fbp_helix_gh(cg, ig, proj, varargin)
%|
%| A single slice rebinning method for cone-beam tomography data
%| collected with a helical source trajectory
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns nt na]	cone-beam projection views (line integrals)
%|
%| option
%|	'short'	1|0		1 for short-scan fan beam (default); 0 for 360
%|	'chat'	1|0		verbosity
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|	sino	[ns na_rebin nz] fan-beam sinogram
%|
%| See helix_example, for information on how to call
%|
%| Equations used are taken from
%| Noo F, Defrise M, Clackdoyle R and Kudo H; Phys. Med. Biol. 44:561-70, 1999
%| "Single-slice rebinning for helical cone-beam CT"
%| @u doi 10.1088/0031-9155/44/2/019
%|
%| Copyright 2010-07-21, Gregory Handy and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test')
	run_mfile_local helix_example
return
end

if nargin < 3, help(mfilename), error(mfilename), end

warn 'this file is obsolete; use rebin_helix and fbp_helix_stat'

arg.short = true;
arg.chat = false;
arg = vararg_pair(arg, varargin);

[img sino] = fbp_helix_do(proj, ig.mask_or, ...
		ig.z, ig.dz, ig.dx, ig.nx, ig.ny, ig.nz, ...
		cg.na, cg.ns, cg.ds, cg.s, cg.t, cg.nt, cg.dt, ...
		cg.dsd, cg.dso, cg.dod, cg.dfs, cg.offset_s, cg.wt, ...
		cg.orbit, cg.orbit_start, cg.pitch, cg.rmax, ...
		cg.source_zs, arg.short, arg.chat);

end % fbp_helix_gh()


%
% fbp_helix_do()
%
function [img sino] = fbp_helix_do(proj, mask2, ...
		zslice, dz, dx, nx, ny, nz, ...
		na, ns, ds, spoints, tpoints, nt, dt, ...
		dsd, dso, dod, dfs, offset_s, wt, ...
		orbit, orbit_start, pitch, rmax, source_zs, short, chat)


% step 1: z-sampling
na1 = na/orbit*360; % # of views in one turn - should be integer!
if abs(round(na1) - na1) > 1e-4
	keyboard
	fail 'bad na/orbit'
end
betas = mod(deg2rad(orbit_start + 360 * [0:(na1-1)]'/na1), 2*pi);

% check to see if the data is used "fully"
if dz > (dso/dsd * dt) * 1.1
	warn('Full use of CB data is not achieved')
end

myPitch = pitch * nt * dso / dsd * dt;
delta = asin(rmax/dso); % fan angle (one sided)

% dist = d in the paper, and allows for either a short scan,
% or for a full 360 scan for each z-slice
if short == 1
	dist = 0.5 * myPitch * (pi + 2 * delta)/(2*pi);
else
	dist = 0.5 * myPitch;

	max_t = (max(abs(spoints)) + dsd^2) / (dso*dsd) * dist;
	if max_t >= max(abs(tpoints))
		error('CT geometry does not allow for 360 degree rebinning')
	end
end


% step 2: rebin the cone-beam data into fan-beam projections
fan_beam_proj = zeros(ns, na1, nz);

% loop over the different view angles
if chat
	printm('Rebinning step beginning')
end
for ia=0:na-1
	ticker(mfilename, ia, na)
	currentZ = source_zs(ia+1);

	upper_limit = currentZ + dist;
	lower_limit = currentZ - dist;

	% acceptable range of z-slices
	iz_list = find((lower_limit <= zslice) & (zslice <= upper_limit));
	if isempty(iz_list), continue, end

	% loop over the acceptable z-slices for the current view angle
	for iz=iz_list'
		% Calculate the values of t to be used for each s
		deltaZ = zslice(iz) - currentZ;
		tpoints = (((spoints).^2 + dsd^2) / (dso*dsd)) * deltaZ; % [ns]

		t_index = tpoints / dt + wt; % [ns]
		it0 = floor(t_index); % [ns]
		it1 = 1 + it0;
		alpha = t_index - it0; % [ns] for linear interpolation
		it0 = it0 + 1; % matlab indexing
		it1 = it1 + 1; % matlab indexing

		if any(it0 < 1) || any(it0 > nt), fail 'bug', end
		if any(it1 < 1) || any(it1 > nt), fail 'bug', end

		% scaling factor due to different ray lengths
		scale = sqrt(spoints.^2 + dsd^2) ...
			./ sqrt(spoints.^2 + tpoints.^2 + dsd^2); % [ns]

		tmp = sub2ind([ns nt], 1:ns, it0') + ia * ns * nt;
		y0 = proj(tmp');
		tmp = sub2ind([ns nt], 1:ns, it1') + ia * ns * nt;
		y1 = proj(tmp');
		fan_beam_proj(:, mod(ia, na1) + 1, iz) = ...
			scale .* (alpha.*y1 + (1-alpha).*y0);
	end
end

if short
	if chat, printm('Preparing sinogram for FBP2'), end

	fb_orbit_start = zeros(nz,1);
	% beginning orbit for each z-slice
	new_orbit = deg2rad(orbit_start);
	zloc = source_zs(1);
	change_betas = betas(2)-betas(1);
	for iz = 1:nz
		orbitFound = 0;
		while orbitFound == 0
			if zslice(iz) < zloc + dist
				fb_orbit_start(iz) = new_orbit;
				orbitFound = 1;
			else
				new_orbit = new_orbit + change_betas;
				zloc = zloc + (source_zs(2) - source_zs(1));
			end
		end
	end

	% # of angles for each z-slice sinogram
%	newNA = floor(2 * dist / myPitch * na1); % GH
	newNA = ceil(2 * dist / myPitch * na1); % JF
	new_fan_beam_proj = zeros(ns, newNA, nz);

	% delete the empty rows of the sinograms
	for iz = 1:nz
		angle = round((fb_orbit_start(iz)-deg2rad(orbit_start))...
				/ change_betas);
		for ia = 0:newNA-1
			nextAngle = angle + ia;
			new_fan_beam_proj(:, ia+1, iz) = ...
			fan_beam_proj(:, mod(nextAngle, na1)+1, iz);
		end
	end

	% same image geometry as before, except for the z direction
	ig = image_geom('nx', nx, 'ny', ny, 'dx', dx, 'mask', mask2);

	img = nans(nx, ny, nz);

	if chat, printm('Performing the filter back projection.'), end
	% filter back project for each z-slice
	for iz = 1:nz
		sg = sino_geom('fan', 'ns', ns, 'na', newNA, ...
			'ds', ds, 'offset_s', offset_s, ...
			'dsd', dsd, 'dod', dod, 'dfs', dfs, ...
			'orbit', 'short', ...
			'orbit_start', rad2deg(fb_orbit_start(iz)));

		if iz == 1 % trick: Parker weights do not depend on orbit_start
			[parker_wt scale180] = fbp_fan_short_wt(sg);
		end

		% apply Parker weighting and scaling
		tmp = scale180 * parker_wt .* new_fan_beam_proj(:,:,iz);
		geomsino = fbp2(sg, ig);
		img(:,:,iz) = fbp2(tmp, geomsino);
	end

	if nargout >= 2
		sino = new_fan_beam_proj;
	end

else % 360
	ig = image_geom('nx', nx, 'ny', ny, 'dx', dx);
	mask2 = true([nx ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1]);
	clear mask2

	img = zeros(nx,ny,nz);

	if chat, printm('Performing the filter back projection.'), end
	% filter back project for each z-slice
	for iz = 1:nz
		sg = sino_geom('fan', 'ns', ns, 'na', na1, ...
				'ds', ds, 'orbit', 360, ...
				'orbit_start', orbit_start, ...
				'offset_s', offset_s, ...
				'dsd', dsd, 'dod', dod, 'dfs', dfs);
		geomsino = fbp2(sg, ig);
		img(:,:,iz) = fbp2(fan_beam_proj(:,:,iz), geomsino);
	end

	if nargout >= 2
		sino = fan_beam_proj;
	end
end

end % fbp_helix_do()
