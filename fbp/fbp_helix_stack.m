  function [img sino] = fbp_helix_stack(cg, ig, sino, orbits, varargin)
%|function [img sino] = fbp_helix_stack(cg, ig, sino, orbits, varargin)
%|
%| Perform slice-by-slice fan-beam FBP.
%| Designed to follow rebin_helix.
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	sino	[ns na nz]	fan-beam sinograms for each slice
%|	orbits	[nz 2]		[orbit orbit_start] for each sinogram
%|
%| option
%|	'window'	''	see fbp2
%|	'short'		1|0	1 for short-scan fan beam (default); 0 for 360
%|	'chat'		1|0	verbosity
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|	sino	[ns na nz]	after applying Parker weighting
%|
%| See helix_example, for information on how to call
%|
%| Copyright 2010-07-21, Gregory Handy and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test')
	run_mfile_local helix_example
return
end

if nargin < 4, help(mfilename), error(mfilename), end

arg.window = ''; % ramp
arg.short = true;
arg.chat = false;
arg = vararg_pair(arg, varargin);

nz = ig.nz;
% 2d version of image geometry
ig = image_geom('nx', ig.nx, 'dx', ig.dx, 'offset_x', ig.offset_x, ...
		'ny', ig.ny, 'dy', ig.dy, 'offset_y', ig.offset_y, ...
		'mask', ig.mask_or);

if any(orbits(:,1) ~= orbits(1,1))
	fail 'only constant orbit done'
end

img = ig.zeros;
sg = sino_geom('fan', 'dsd', cg.dsd, 'dso', cg.dso, ...
	'ns', cg.ns, 'ds', cg.ds, 'offset_s', cg.offset_s, ...
	'orbit', orbits(1,1), ...
	'na', size(sino, 2)); % # views in fan-beam sinogram, not helix cg.na

% trick: Parker weights do not depend on orbit_start
if arg.short
	[parker_wt scale180] = fbp_fan_short_wt(sg);
else % 360
	parker_wt = 1;
	scale180 = 1;
end

for iz = 1:nz % each slice
	ticker(mfilename, iz, nz)

	sg.orbit_start = orbits(iz,2); % sino_geom for this slice!

	% apply Parker weighting and scaling
	sino(:,:,iz) = scale180 * parker_wt .* sino(:,:,iz);
	geom = fbp2(sg, ig);
	img(:,:,iz) = fbp2(sino(:,:,iz), geom, 'window', arg.window);
end

end % fbp_helix_stack()
