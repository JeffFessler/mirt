% test_denoise_quad.m
% Test denoising with quadratic regularization.
% motivated by sense map issues with michael allison
% Copyright 2005, Jeff Fessler, University of Michigan
% 2011-06-14: mex version and mat version differ due to cgrad,diff edge issues

if ~isvar('R'), printm 'R'
	f.l2b = -2;
	f.order = 2;

	ig = image_geom('nx', 2^5, 'ny', 2^4-2, 'dx', 1);
	if 0 % experiment with excluding image boundary
		ig.mask(:,[1 end]) = false;
		ig.mask([1 end],:) = false;
	end

	args = {ig.mask, 'beta', 2^f.l2b, 'order', f.order};
	R = Reg1(args{:});
%	R = Reg1(args{:}, 'type_penal', 'mat', 'type_diff', 'spmat');
%	R = Reg1(args{:}, 'type_penal', 'mat'); % this one is fine
	psf = qpwls_psf(1., R, 1., ig.mask); % expected blur (at image center)
end

f.niter = 100;
% todo: changes a lot with iteration so make a movie to watch it!

if 1 % run qpwls algorithm for regularized fitting
	A = 1;
	init = ig.zeros;
	yy = ig.unitv('c', [ig.nx/2-1 1]);
	xq = qpwls_pcg1(init(ig.mask), A, 1, yy(:), R.C, ...
		'stop_grad_tol', 1e-8, ...
		'niter', f.niter, 'isave', 'last', 'chat', 1);
	xq = ig.embed(xq);

	xp = pwls_pcg1(init(ig.mask), A, 1, yy(:), R, ...
		'stop_grad_tol', 2e-7, ...
		'niter', f.niter, 'isave', 'last', 'chat', 1);

	xp = ig.embed(xp);
	im plc 2 2
	im(1, psf)
	im(2, xq)
	im(3, xp)
	im(4, xp-xq)
	equivs(xq,xp) % fails for mex version, due to boundary issues
return
end

if 1 % qpwls algorithm for case with 0 weights outside ROI
	A = 1;
	tmp = abs(ig.xg) < ig.nx/3 & abs(ig.yg) < ig.ny/4; % flat
%	tmp = tmp & (abs(ig.xg) > ig.nx/4 & abs(ig.yg) > ig.ny/5); % box
	tmp = (1 + cumsum(tmp, 1) / ig.nx) .* tmp; % ramp
	tmp = tmp .* exp(1i * tmp * 0.5); % make complex!
	yy = tmp;
	ww = abs(tmp) > 0; % weighting
	im pl 2 3
	im clf, im(1, tmp), cbar
	im(2, ww), cbar
	im(3, psf), cbar
	drawnow

	W = diag_sp(ww(:));
	init = ig.zeros;
	xh = qpwls_pcg1(init(ig.mask), A, W, ...
		yy(:), R.C, 'niter', f.niter);
	xh = ig.embed(xh);
	im(4, xh), cbar

	im subplot 5
	plot(	...
		ig.x, real(tmp(:,end/2)), 'r-', ...
		ig.x, real(xh(:,end/2)), 'r.', ...
		ig.x, imag(tmp(:,end/2)), 'g-', ...
		ig.x, imag(xh(:,end/2)), 'g.')
end



return % not done below here ............

if ~isvar('C'), printm 'C'
	C = R.C;
	Cf = C(:,:);
	Ch = Cf' * Cf;
end
	A = eye(ig.nx);

return

if arg.slow, printm 'make sparse hessian'
	A = spdiag(arg.bodycoil(mask(:)), 'nowarn');
	C = R.C;
	H = A' * A + C' * C;
	if 0 % test equivalence of gradients
		Ro = Reg1(args{:});
		t1 = C' * (C *  sinit(:));
		t2 = Ro.cgrad(Ro, sinit(:));
		equivs(t1, t2)
	end
else
	A = diag_sp(arg.bodycoil(mask(:)));
end

if streq(arg.isave, 'all')
	smap = zeros(nx,ny,ncoil,arg.niter+1);
else
	smap = zeros(nx,ny,ncoil,length(arg.isave));
end

if streq(arg.precon, 'diag')
	if R.order == 2
		arg.precon = 6 * sum(R.beta); % (-1)^2 + 2^2 + (-1)^2
	else
		arg.precon = 2 * sum(R.beta); % (-1)^2 + (1)^2
	end
	arg.precon = abs(arg.bodycoil).^2 + arg.precon;
%	im(arg.precon), prompt
	arg.precon = 1 ./ arg.precon;
	arg.precon = diag_sp(arg.precon);
end

for ic = 1:ncoil
	ytmp = ykj(:,:,ic);
	if arg.slow % numerical inverse for testing
		ytmp = double(ytmp);
		tmp = A' * ytmp(:);
		tmp = H \ tmp;

		if nargout > 2 % cost
			cost(:,ic) = pwls_cost(tmp, A, 1, ytmp(mask(:)), R);
		end
	end
	smap(:,:,ic,:) = embed(tmp, mask);
end
