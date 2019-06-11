%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 4-level multiresulution nonrigid image registration with Jacobian penalty
%% Conjugate gradient with Newton's first step was used
%% Cubic B-spline bases are used for both image and deformation models
%%
%% Copyright July 2007, Se Young Chun and Jeff Fessler, University of Michigan
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear;

if ~isvar('Xctsrc')
	% generate a source image, a circle
	Xctsrc = zeros(256);
	for i = 1 : 256
		for j = 1 : 256
			Xctsrc(i, j) = ((i-128)^2 + (j-128)^2 < 60^2);
		end
	end
	Xctsrc = single(Xctsrc); % All data should be single type

	% generating a target image, a rectangle
	Xcttar = zeros(256);
	Xcttar([128-60:128+60], [128-60:128+60]) = ones(121);
	Xcttar = single(Xcttar);
end


if ~isvar('kg')
	% image information, source and target images have same info
	ig = image_geom('nx', 256, 'dx', 1, 'ny', 256, 'dy', 1);
	% knot (deformation) information
	kg = knot_geom('nx', 32,  'mx', 8, 'offset_x', 0, ...
			'ny', 32, 'my', 8, 'offset_y', 0);
end


if ~isvar('alphax')
	% B spline coefficients for source image
	Coeffsrc = BsplVal2CoMirr(Xctsrc);

	% generating deformation matrices B using image and knot information
	Bx = makeB(ig, kg);
	By = Bx;

	% initializing deformation coefficients, should be of single type
	alphax = single(zeros(kg.dim));
	alphay = alphax;
end


if ~isvar('beta')
	% First conjugate gradient iteration
	[W Wgx Wgy] = makeW({Bx, By}, {alphax, alphay});

	diff = W*Coeffsrc - Xcttar(:);
	iter = 1;
	nd(iter) = norm(diff);
	printf('Iter%d, norm(diff) = %d', iter, nd(iter));

	gsx = Wgx*Coeffsrc;
	gsy = Wgy*Coeffsrc;

	% gradients
	gx = (Bx'*(gsx.*diff));
	gy = (By'*(gsy.*diff));

	% conjugate directions
	dx = -gx;
	dy = -gy;

	tx = Bx*dx;
	ty = By*dy;

	% Newton method first iteration step size
	gamma = - diff'*(gsx.*tx + gsy.*ty) / sum( (gsx.*tx + gsy.*ty).^2 );

	%deformation update
	alphax = alphax + gamma * reshape(dx, kg.dim);
	alphay = alphay + gamma * reshape(dy, kg.dim);

	% Conjugate gradient iterations (max 300 iterations)
	while ((norm(diff) > 1)&&(iter < 300))
		iter = iter + 1;
		[W Wgx Wgy] = makeW({Bx, By}, {alphax, alphay});

		diff = W*Coeffsrc - Xcttar(:);
		nd(iter) = norm(diff);
		printf('Iter%d, norm(diff) = %d', iter, nd(iter));

		gx_old = gx;
		gy_old = gy;

		gsx = Wgx*Coeffsrc;
		gsy = Wgy*Coeffsrc;

		gx = (Bx'*(gsx.*diff));
		gy = (By'*(gsy.*diff));

		beta = ( sum(gx.*(gx-gx_old)) + sum(gy.*(gy-gy_old)) ) / ...
			( sum(gx_old.*gx_old) + sum(gy_old.*gy_old) );
		beta = max(beta, 0);

		dx = -gx + beta * dx;
		dy = -gy + beta * dy;

		tx = Bx*dx;
		ty = By*dy;

		gamma = - diff'*(gsx.*tx + gsy.*ty) / sum( (gsx.*tx + gsy.*ty).^2 );

		alphax = alphax + gamma * reshape(dx, kg.dim);
		alphay = alphay + gamma * reshape(dy, kg.dim);
	end
end


% Display results
im plc 2 2
im(1, Xctsrc), cbar
title 'source image'

[W Wgx Wgy] = makeW({Bx, By}, {alphax, alphay});
y = reshape(W*Coeffsrc, ig.dim);
im(2, y), cbar
title 'deformed image'

im(3, Xcttar - y), cbar
title 'diff image'

im(4, W.handle_ufun(W));
title 'deformation grid'
