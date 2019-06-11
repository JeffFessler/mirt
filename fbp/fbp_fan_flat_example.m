% fbp_fan_flat_example
% example of how to use fbp_fan_flat.m
% Copyright Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'specify geometry and synthesize sinogram'
	sg = sino_geom('fan', 'nb', 280, 'na', 400, 'ds', 4, ...
		'orbit', 360, ...
		'offset_s', 0.25, 'dsd', 949, 'dod', 408, 'dfs', inf); 

	ell = [20 -15 200 200 30 10; 80 0 9 9 0 9; -70 -40 40 40 0 -8];
	sino = ellipse_sino(sg, ell, 'oversample', 3);

	ig = image_geom('nx', 256, 'ny', 252, 'dx', 2, 'down', 2);
	xtrue = ellipse_im(ig, ell, 'oversample', 4);

	im plc 2 3
	im(1, xtrue, 'x'), cbar
	im(2, sino, 'sino'), cbar
prompt
end

% fan-beam reconstruction
if ~isvar('recon'), printm 'fbp recon'
	geom = fbp2(sg, ig);
	recon = fbp2(sino, geom);
end

if 0 % verify consistency with old way
	G = Gtomo2_dscmex(sg, ig); % system object
	G = struct(G);
	recon0 = fbp_fan_flat(sino, G, 'ramp');
	mask = recon0 ~= 0;
	recon = recon .* mask;
	max_percent_diff(recon, recon0)
	im clf, im([recon recon0])
return
end

if im
	im(4, recon, 'recon'), cbar
	im(5, recon - xtrue, 'error'), cbar
	iy = ig.ny/2; ix = 1:ig.nx;
	subplot(133)
	plot(ix, xtrue(ix,iy), '-', ix, recon(ix,iy), '--')
	axis([1 ig.nx -0.5 20])
	legend('true', 'recon', 'location', 'south')
end
