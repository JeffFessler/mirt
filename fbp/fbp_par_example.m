% fbp_par_example.m
% FBP parallel-ray example using fbp_dsc.m
% Copyright 2005-8-10, Jeff Fessler, The University of Michigan

redo = 1;
if redo || ~isvar('sino')
	down = 2;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	sg = sino_geom('par', 'nb', 888, 'na', 984, 'dr', 541/949, ...
		'orbit', 180, ... % curiously, 360 has less aliasing artifact
		'offset_r', 0.25, 'orbit_start', 0, 'down', down);

	ell = [20 10 200 200 0 1; 38 27 10 10 0 1];
%	ell = [0 0 200 200 0 0; 13 13 10 10 0 1000];

	sino = ellipse_sino(sg, ell, 'oversample', 4);
	xtrue = ellipse_im(ig, ell, 'oversample', 2);

	im pl 2 3
	im(1, xtrue, 'x'), cbar
	im(2, sino, 'sino'), cbar
prompt
end

% fbp reconstruction
if redo || ~isvar('recon')
	if 1
		args = aspire_pair(sg, ig, 'system', 9);
		recon = fbp_dsc(sino, 1, args);
	else
		tmp = fbp2(sg,ig);
		recon = fbp2(sino, tmp);
	end
end

if 0 % to examine aliasing
	clim = [9.5 10.5];
	clim = [-200 200];
	clim = [-50 50];
	clim = [0.8 1.2];
	clim = [0.9 1.1];
	im(4, recon, clim), cbar
%	im(recon .* (xtrue == 0)), cbar
return
end
%im(4, recon, 'recon'), cbar
im(5, recon - xtrue, 'error'), cbar
iy = ig.ny/2; ix = 1:ig.nx;
subplot(133)
plot(ix, xtrue(ix,iy), '-', ix, recon(ix,iy), '--')
axis([1 ig.nx -0.5 1.1 * max(xtrue(:))])
legend('true', 'recon', 'location', 'south')
