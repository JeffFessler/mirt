% pwls_sps_zoom.m
% An example of using PWLS-SQS algorithm for edge-preserving image zooming.
% Copyright Jan 2001, Jeff Fessler, University of Michigan

% generate data
if ~isvar('yi')
	ig = image_geom('nx', 110, 'ny', 128, 'dx', 1);
	xtrue = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	im plc 2 2
	im(1, xtrue, 'xtrue (high res)')

	% generate sparse "down sampling" system matrix
	down = 2;
	n.b = ig.nx / down;
	n.a = ig.ny / down;
	if 1
		A = Gdown([ig.nx ig.ny], 'down', down);
	else
		Ax = kron(speye(n.b), sparse(ones(1,down)/down));
		Ay = kron(speye(n.a), sparse(ones(1,down)/down));
		A = kron(Ay, Ax); clear Ax Ay
		A = Gsparse(A, 'mask', ig.mask, 'odim', [n.b n.a]);
	end

	% y = A * x does the down sampling in matrix form!
	yi = A * xtrue;
	im(2, yi, 'yi (low res)')
prompt
end


% elementary 0-order hold zooming and cubic interpolation
if ~isvar('xc')
	x0 = down^2 * ig.shape(A' * yi(:));
	im(3, x0, 'x0, 0-order hold')

	xc = interp2(yi', 0.5+[0:(n.b*down-1)]'/down, 0.5+[0:(n.a*down-1)]/down)';
	xc(isnan(xc)) = 0;
	im(4, xc, 'xc, cubic interpolation')
prompt
end


% regularization
if ~isvar('Rn')
	Rq = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', 'quad', 'beta', 2^(-4));
	Rn = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', 'huber', 'beta', 2^(-7), 'delta', 0.3);
prompt
end


% SQS with quadratic penalty
if ~isvar('xq')
	f.niter = 40;
	xq = pwls_sqs_os(xc(:), A, yi, Rq, 'niter', f.niter);
	xq = ig.embed(xq);

	im clf, im(xq, 'QPWLS-SQS iterations')
prompt
end


% SQS with nonquadratic edge-preserving penalty
if ~isvar('xs'), printm 'xs'
	f.niter = 50;
	xs = pwls_sqs_os(xq(:), A, yi, Rn, 'niter', f.niter-1, 'isave', 'all');
	xs = ig.embed(xs);

	im clf, im(xs, 'PWLS-SQS iterations')
prompt
end

if im
	rms = @(x) sqrt(mean(abs(x).^2));
	r0 = rms(col(xtrue-x0));
	rc = rms(col(xtrue-xc));
	rq = rms(col(xtrue-xq));
	rs = rms(col(xtrue-xs(:,:,end)));

	im plc 2 3
	clim = [0 6];
	im(1, xtrue, clim, 'truth'), cbar
	im(2, yi, clim, 'yi'), cbar
	axis([1 ig.nx 1 ig.ny])
	im(3, x0, clim, 'x0, zero-order'), cbar
	xlabel(num2str(r0))
	im(4, xc, clim, 'xc, cubic'), cbar
	xlabel(num2str(rc))
	im(5, xq, clim, 'SQS quadratic'), cbar
	xlabel(num2str(rq))
	im(6, xs(:,:,end), clim, 'SQS edge preserving'), cbar
	xlabel(num2str(rs))
end
