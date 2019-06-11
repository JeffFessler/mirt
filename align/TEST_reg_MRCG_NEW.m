%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                            
%% 4-level multiresulution nonrigid image registration with New penalty
%% Conjugate gradient with Newton's first step was used 
%% Cubic B-spline bases are used for both image and deformation models   
%%                                                                            
%% Copyright April 2008, Se Young Chun and Jeff Fessler, University of Michigan
%%                                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isvar('alphay'), printm 'registering images'
	betaJ = 100;

	%kg.mx will be set 8
	mx = 0.99*4*0.5; % detJ = 0.01
	Mx = 9*4; % detJ = 100.2450
	% New penalty functions
	pm0 = @(x) (x+mx).^2.*(x<=-mx)/2 + (x-mx).^2.*(x>mx)/2;
	pm1 = @(x) (x+mx).*(x<=-mx) + (x-mx).*(x>mx);
	pm2 = @(x) (x<=-mx) + (x>mx);
	pM0 = @(x) (x+mx).^2.*(x<=-mx)/2 + (x-Mx).^2.*(x>Mx)/2;
	pM1 = @(x) (x+mx).*(x<=-mx) + (x-Mx).*(x>Mx);
	pM2 = @(x) (x<=-mx) + (x>Mx);

	% generating a source image, a circle
	Xctsrc = zeros(256);
	for i = 1 : 256
		for j = 1 : 256
			Xctsrc(i, j) = ((i-128)^2 + (j-128)^2 < 60^2);
		end
	end
	Xctsrc = single(Xctsrc);

	% generating a target image, a rectangle
	Xcttar = zeros(256);
	Xcttar([128-60:128+60], [128-60:128+60]) = ones(121);
	Xcttar = single(Xcttar);

	% generating image pyramid for multiresolution image registration
	Xdnsrc2 = BsplReduce(Xctsrc);
	Xdnsrc4 = BsplReduce(Xdnsrc2);
	Xdntar2 = BsplReduce(Xcttar);
	Xdntar4 = BsplReduce(Xdntar2);
	Xdnsrc8 = BsplReduce(Xdnsrc4);
	Xdntar8 = BsplReduce(Xdntar4);

	iter = 0;
	for rate = [8 4 2 1]

		% image information, source and target images have same info
		ig = image_geom('nx',256/rate,'dx',1,'ny',256/rate,'dy',1);
		% knot (deformation) information
		kg = knot_geom('nx', 64/rate,  'mx', 4, 'offset_x', 0, ...
				'ny', 64/rate, 'my', 4, 'offset_y', 0);

		% Difference matrices for New penalty
        	Px = CsparsePx(kg.dim);
        	Py = CsparsePy(kg.dim);

		% generating deformation matrices B using image and knot info
		B = makeB(ig, kg);

		iter = iter + 1;

		switch(rate)
			case 8
				% smallest images
				Xsrc = Xdnsrc8;
				Xtar = Xdntar8;

				% B spline coefficients for source image
				Csrc = BsplVal2CoMirr(Xsrc);

				% initializing deformation coefficients, 
				% should be of single type
				alphax = single(zeros(kg.dim));
				alphay = alphax;
			case 4
				% medium size images
				Xsrc = Xdnsrc4;
				Xtar = Xdntar4;

				% B spline coefficients for source image
				Csrc = BsplVal2CoMirr(Xsrc);

				% converting deformation info to the higher 
				% resolution deformation info
				tmp = BsplCo2ValZero(alphax, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphax = BsplVal2CoZero(tmp);
				tmp = BsplCo2ValZero(alphay, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphay = BsplVal2CoZero(tmp);
			case 2
				% medium size images
				Xsrc = Xdnsrc2;
				Xtar = Xdntar2;

				% B spline coefficients for source image
				Csrc = BsplVal2CoMirr(Xsrc);

				% converting deformation info to the higher 
				% resolution deformation info
				tmp = BsplCo2ValZero(alphax, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphax = BsplVal2CoZero(tmp);
				tmp = BsplCo2ValZero(alphay, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphay = BsplVal2CoZero(tmp);
			case 1
				% largest images
				Xsrc = Xctsrc;
				Xtar = Xcttar;

				% B spline coefficients for source image
				Csrc = BsplVal2CoMirr(Xsrc);

				% converting deformation info to the higher 
				% resolution deformation info
				tmp = BsplCo2ValZero(alphax, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphax = BsplVal2CoZero(tmp);
				tmp = BsplCo2ValZero(alphay, kg.dim, ...
					[0 0], [2 2], {});
				tmp = 2*tmp;
				alphay = BsplVal2CoZero(tmp);
		end

		% First conjugate gradient iteration for each resolution
		[W Wgx Wgy] = makeW({B, B}, {alphax, alphay});

		diff = W*Csrc - Xtar(:);
		nd(iter) = rate*norm(diff);
		printf('==== Iter%d, norm(diff) = %d', iter, nd(iter));

		gsx = Wgx*Csrc; gsy = Wgy*Csrc;

		pxx = Px*double(alphax(:)); pyx = Py*double(alphax(:));
        	pxy = Px*double(alphay(:)); pyy = Py*double(alphay(:));

		% gradients
		gx = (B'*(gsx.*diff)) + betaJ*single(Px'*pM1(pxx)+Py'*pm1(pyx));
		gy = (B'*(gsy.*diff)) + betaJ*single(Px'*pm1(pxy)+Py'*pM1(pyy));

		% conjugate directions
		dx = -gx; dy = -gy;

		tx = B*dx; ty = B*dy;

        	pxdx = Px*double(dx(:)); pydx = Py*double(dx(:));
        	pxdy = Px*double(dy(:)); pydy = Py*double(dy(:));

		% Newton method first iteration step size
		gamma_num = gx'*dx + gy'*dy;
		gamma_den = sum( (gsx.*tx + gsy.*ty).^2 ) + betaJ*( ...
			(pxdx.^2)'*(pM2(pxx)) + (pydx.^2)'*(pm2(pyx)) ...
                	+ (pxdy.^2)'*(pm2(pxy)) + (pydy.^2)'*(pM2(pyy)) );
		gamma = - gamma_num / gamma_den;

		%deformation update	
		alphax = alphax + gamma * reshape(dx, kg.dim);
		alphay = alphay + gamma * reshape(dy, kg.dim);

		% Conjugate gradient iterations (max 100 iterations)
		flag = 0;
        	subiter = 0;
        	while (flag < 0.5)
			iter = iter + 1;
			subiter = subiter + 1;
			[W Wgx Wgy] = makeW({B, B}, {alphax, alphay});

			diff = W*Csrc - Xtar(:);
			printf('Iter%d, norm(diff) = %d',iter,rate*norm(diff));
			nd(iter) = rate*norm(diff);

			gx_old = gx; gy_old = gy;
	
			gsx = Wgx*Csrc; gsy = Wgy*Csrc; 

			pxx = Px*double(alphax(:)); pyx = Py*double(alphax(:));
        		pxy = Px*double(alphay(:)); pyy = Py*double(alphay(:));

			% gradients
			gx = (B'*(gsx.*diff)) + betaJ*single(...
					Px'*pM1(pxx)+Py'*pm1(pyx));
			gy = (B'*(gsy.*diff)) + betaJ*single(...
					Px'*pm1(pxy)+Py'*pM1(pyy));

			grad = max( sqrt( gx.^2 + gy.^2) );

			if ( grad < 10^(-5) ) | (subiter > 100)
                        	flag = 1;
                	end

			beta = (sum(gx.*(gx-gx_old))+sum(gy.*(gy-gy_old)))/...
				(sum(gx_old.*gx_old)+sum(gy_old.*gy_old));
			beta = max(beta, 0);

			dx = -gx + beta * dx;
			dy = -gy + beta * dy;
	
			tx = B*dx;
			ty = B*dy;

        		pxdx = Px*double(dx(:)); pydx = Py*double(dx(:));
        		pxdy = Px*double(dy(:)); pydy = Py*double(dy(:));

			gamma_num = gx'*dx + gy'*dy;
			gamma_den = sum( (gsx.*tx + gsy.*ty).^2 ) + betaJ*( ...
				(pxdx.^2)'*(pM2(pxx)) + (pydx.^2)'*(pm2(pyx))...
                		+ (pxdy.^2)'*(pm2(pxy))+(pydy.^2)'*(pM2(pyy)));
			gamma = - gamma_num / gamma_den;

			alphax = alphax + gamma * reshape(dx, kg.dim);
			alphay = alphay + gamma * reshape(dy, kg.dim);
		end
	end
end

% Display results
im pl 2 3
im(1, Xctsrc), cbar
title 'source image'

[W Wgx Wgy] = makeW({B, B}, {alphax, alphay});
y = reshape(W*Csrc, ig.dim);
im(2, y), cbar
title 'deformed image'

im(3, Xcttar - y), cbar
title 'diff image'

im(4, W.handle_ufun(W))
title 'deformation grid'

[B Bgx Bgy] = makeB(ig, kg);
xx = Bgx*alphax; yy = Bgy*alphay;
xy = Bgx*alphay; yx = Bgy*alphax;
detJ = reshape((1+xx).*(1+yy) - xy.*yx, ig.dim);

im(5, detJ>0)
title 'detJ - binary'

if 1 printm 'checking finer jacobian'
        if ~isvar('detJm')
		igm = image_geom('nx', 2560, 'dx', 1, 'ny', 2560, 'dy', 1);
		kgm = knot_geom('nx', 64,  'mx', 40, 'offset_x', 0, ...
       	         		'ny', 64, 'my', 40, 'offset_y', 0);
		[Bm Bmgx Bmgy] = makeB(igm, kgm);
		xx = 10*(Bmgx*alphax); yy = 10*(Bmgy*alphay);
		xy = 10*(Bmgx*alphay); yx = 10*(Bmgy*alphax);
		detJm = reshape((1+xx).*(1+yy) - xy.*yx, igm.dim);
	end

	im(6, detJm>0)
	title 'detJ - binary (finer)'
end
