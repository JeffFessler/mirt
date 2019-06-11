% ir_fbp2_test_non180.m
% Test FBP with non-180 orbit
% Copyright 2014-12-03, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'sino'
	down = 2;
	ig = image_geom('nx', 488, 'ny', 512, 'fov', 500);
	ig = ig.downsample(down);

	% parallel-beam
	sg = sino_geom('par', 'nb', 888, 'na', 984, ...
		'orbit_start', -400, ... % strest test!
		'orbit', -220, ... % stress test!
		'strip_width', 'd', ...
		'dr', 541/949, 'offset_r', 0.25);
	sg = sg.downsample(down);

	clim = (1 + [-1 1] * 0.05) * 1000;
	[xtrue ell] = ellipse_im(ig, [], 'oversample', 2, 'hu_scale', 1000);
	sino = ellipse_sino(sg, ell, 'oversample', 2);

	im plc 2 3
	im(1, xtrue, 'xtrue', clim), cbar
	im(4, sino, 'sino'), cbar

	t = floor(min(sg.nb * sg.d, ig.fov)/2);
	ig.mask = ellipse_im(ig, [0 0 t t 0 1]) > 0;
	im(6, ig.mask)
%prompt
end

% conventional FBP reconstruction
if ~isvar('r_std'), printm 'fbp std'
	fg.std = fbp2(sg, ig);
	cpu etic
	r_std = fbp2(sino, fg.std);
	cpu etoc 'fbp std recon time'

	im(2, r_std, 'FBP matlab std', clim), cbar
	im(5, r_std - xtrue, 'error', [-1 1]*100), cbar
%prompt
end
