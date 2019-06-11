% ir_ip_unwrap_example.m
% example of phase unwrapping using optimization transfer
% Copyright 2003-6-15, Jeff Fessler, University of Michigan

if ~isvar('x.true')
%	x.max_diff = pi/2;
	x.max_diff = 8*pi/9;
	nx = 41; ny = 31;
	p.x = [1:nx]-(nx/2+1);
	p.y = [1:ny]-(ny/2+1);
	x.true = gaussian_kernel(nx/3, (nx-1)/2) ...
		* gaussian_kernel(nx/3, (ny-1)/2)';
	% create phase with certain maximum
	x.true = x.true * x.max_diff / ...
	max(max(col(abs(diff(x.true,1)))), max(col(abs(diff(x.true,2)))));

	% true magnitude
	rng(0)
	mag = 0.2 + 0.8 * rand(nx, ny); % random
	mag(1:round(nx/2),:) = 1; % uniform
	im plc 2 4
	x.clim = [0 9*pi];
	im(1, x.true, 'x true', x.clim), cbar
	im(2, mag, 'mag'), cbar
end

if ~isvar('yi'), printm 'simulate noisy data'
	f.sig = 0.3;
	rng(0)
	yi = mag .* exp(1i * x.true) + f.sig * randn(size(x.true));
	x.u = unwrap(unwrap(angle(yi))')';
	x.u = 0.5 * (x.u + unwrap(unwrap(angle(yi)')'));

	im(3, abs(yi), '|y|'), cbar
	im(4, angle(yi), '\angle y'), cbar
%	im(4, phase(yi), '\angle y'), cbar
	im(5, x.u, 'matlab unwrap', x.clim), cbar
end

if 1 || ~isvar('x.sqs'), printm 'iterative estimation'
	f.l2b_q = -3; % choosing this beta may be the hardest part!
	Rq = Reg1(true(nx,ny), 'type_denom', 'matlab', 'beta', 2^f.l2b_q);

	niter = 200;
	niter = 40;
	x.init = x.u;
%	x.init = angle(yi);
	x.sqs = ir_ip_unwrap_sqs(x.init(:), yi(:), Rq, niter);
	x.sqs_end = reshape(x.sqs(:,end), nx, ny);

	im(6, x.sqs_end, 'fessler sqs unwrap', x.clim), cbar
end

cost = abs(yi).^2 .* (1 - cos(x.sqs_end - angle(yi)));
im(7, cost, 'data mismatch'), cbar

a = angle(yi(:));
cost = abs(yi(:)).^2' * (1 - cos(x.sqs - a(:,ones(1,niter))));
cost = cost + Rq.penal(Rq, x.sqs).';
if im
	im subplot 8, plot(0:niter-1, cost, '-o')
	xlabel 'iteration', title 'cost'
end
