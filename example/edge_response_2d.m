% edge_response_2d.m
% examine variation in edge response for 2D CT
% Copyright 2011-02-23, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'setup geometry, image, sinogram'
	down = 2;
	ig = image_geom('nx', 512, 'fov', 50, 'down', down);
	ig.mask = ig.circ > 0;
	sg = sino_geom('ge1', 'units', 'cm', 'down', down);

	% system object
	A = Gtomo2_dscmex(sg, ig);

	% image
	
	steps = [1 2 4 8] * 10;
        nstep = numel(steps);
	ell = [0 0 15 10 0 1000];
	leg = {};
	for ii=1:nstep
		ell = [ell; [-15+ii*6 0 2 2 0 steps(ii)]];
		leg{end+1} = sprintf('%d HU', steps(ii));
	end
	xtrue = ellipse_im(ig, ell, 'oversample', 2);

	im plc 2 2
	clim = [800 1200];
	im(1, xtrue, 'x', clim), cbar

%	sino = ellipse_sino(sg, ell, 'oversample', 2);
	sino = A * xtrue; % cheat
	im(2, sino, 'sino'), cbar
prompt
end


if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, ig);
	fbp = fbp2(sino, tmp, 'window', 'hanning,0.5');
	im(3, fbp, 'FBP', clim), cbar
prompt
end


if ~isvar('kappa'), printm 'kappa'
	wi = ones(size(sino)); % no statistical weighting
	kappa = ig.mask;
end


% use local psf to help select beta
if ~isvar('R'), printm 'R'
	f.l2b = 4; % maybe a bit too big, but ok for now
	f.delta = 10;
	R = Reg1(kappa, 'beta', 2^f.l2b, 'pot_arg', {'hyper3', f.delta});
	qpwls_psf(A, R, 1, ig.mask, Gdiag(wi), 'loop', 1);
prompt
end


if ~isvar('xpwls'), printm 'iterative reconstruction'
	Ab = Gblock(A, 41); % 41 subsets
	f.niter = 20;
	xpwls = pwls_sps_os(fbp(ig.mask), sino, wi, Ab, R, f.niter);
	xpwls = ig.embed(xpwls(:,end));
	im(4, xpwls, clim), cbar
prompt
end

im plc 3 3
pl = @(i) subplot(330+i);

im(1, 'notick', ig.x, ig.y, fbp, clim, 'FBP')
im(4, 'notick', ig.x, ig.y, xpwls, clim, 'PWLS')

tmp1 = fbp(:,end/2+1);
tmp2 = xpwls(:,end/2+1);

pl(2)
plot(ig.x, tmp1)
axis([-20 20 0 1200])
title 'Profile'
xtick([-1 0 1] * 15), ytick([0 1000])

pl(5)
plot(ig.x, tmp2)
axis([-20 20 0 1200])
title 'Profile'
xtick([-1 0 1] * 15), ytick([0 1000])

xx = outer_sum(ig.x, -ell(2:end,1));
pl(3)
plot(xx, (tmp1-1000) * (1 ./ steps), '.-')
ax = [-3 -1 -0.25 1.25];
axis(ax), legend(leg{:}, 'location', 'seo')
title 'Edge response'
xtick([-3 -2 -1]), ytick([0 1])

pl(6)
plot(xx, (tmp2-1000) * (1 ./ steps), '.-')
axis(ax), legend(leg{:}, 'location', 'seo')
title 'Edge response'
xtick([-3 -2 -1]), ytick([0 1])

% ir_savefig eps_c edge_response_2d
