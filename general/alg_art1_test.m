% alg_art1_test.m
% test ART etc.
% Copyright 2006-4-2, Jeff Fessler, University of Michigan

% generate data
if ~isvar('yi'), printm 'yi'
	ig = image_geom('nx', 16, 'ny', 14, 'dx', 1);
	ig.mask = ig.circ > 0;
	sg = sino_geom('par', 'nb', ig.nx, 'na', round(0.6 * ig.nx / 2)*2, ...
		'dr', 1);
	A = Gtomo2_strip(sg, ig);
	xtrue = ig.circ(5);
	yi = A * xtrue;
end


% ART1
if ~isvar('xart1'), printm 'art1'
	f.niter = 34;
	xinit = ig.zeros;

	xart1 = alg_art1(xinit(ig.mask), A', yi, 'niter', f.niter, 'isave', 'all');
	xart1 = ig.embed(xart1);
	im clf, im(xart1, 'ART1 iterates')
prompt
end

if ~isvar('yp')
	yp = A * xart1;
	resid = repmat(yi, [1 1 f.niter+1]) - yp; 
end

	tmp = reshape(resid, [], f.niter+1);

if im
	im plc 3 3

	im(1, xtrue, 'xtrue')
	im(2, xart1, 'ART1')
	im(3, ig.mask, 'mask')
	im(4, yi, 'yi'), cbar
	im(5, yp(:,:,end), 'yp'), cbar
	im(6, yp(:,:,end)-yi, 'yp-yi'), cbar
	im(7, xart1(:,:,end), 'last iter'), cbar
	im(8, xart1(:,:,end)-xtrue, 'err'), cbar
	im subplot 9
	semilogy(0:f.niter, sqrt(mean(tmp.^2)))
end
