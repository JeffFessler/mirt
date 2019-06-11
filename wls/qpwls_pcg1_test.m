% qpwls_pcg1_test
% Test qpwls_pcg1_test with tiny (complex) example having exact solution.
%
% Copyright 2005-9-14, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), printm 'data'
	nx = 20;
	ny = 18;
	mask = true(nx,ny);
	%mask(nx) = false;
	rng(0)
	np = sum(mask(:));
	G = sparse(find(mask(:)), 1:np, ones(np,1), nx*ny, np);
	gg = 7 + randn(nx,ny) + 2i * randn(nx,ny);
	gg(end-4:end,:) = 0;
	gg(end/2,:) = 0;
	gg(1:3,:) = 0;
	G = sparse(1:nx*ny, 1:nx*ny, gg(:), nx*ny, nx*ny) * G;
	G = full(G);
	if 0
		G = sparse([1:nx*ny]+cumsum(~mask), 1:np, gg);
		G = diag(1 + 0.2 * randn(nx*ny,1)); G(~mask(:),:);
	end
	W = 1/mean(abs(gg(:)))^2;
	x = [1:nx]'*[1:ny] + nx*0.5i;
	yb = reshape(G * x(mask), nx, ny);
	yi = yb + randn(size(yb)) + 1i * randn(size(yb));
	im clf, pl = 320;
	im(pl+1, abs(x), 'x true'), cbar
	im(pl+2, abs(yi), 'yi'), cbar
	im(pl+5, abs(gg), 'gg'), cbar
prompt
end


if ~isvar('R'), printm 'make R'
	R = Robject(mask, 'beta', 2^0, 'order', 2, 'offsets', 1);
	C = R.C;
	C = C(:,:);
	im(pl+3, C', 'C')

	qpwls_psf(G, R.C, 1, mask, W); % check resolution, adjust beta if needed
prompt
end


%
% exact solution
%
if ~isvar('xhat'), printm 'exact xhat'
	F = full(G' * W * G);
%	C1 = R.C1(:,:);
%	Rh = full(C1' * spdiag(R.wt) * C1);	% penalty Hessian
	Rh = C' * C;	% penalty Hessian
	printm('cond F = %g', cond(F))
	printm('cond(F+R) = %g', cond(F+Rh))
	xhat = (F + Rh) \ (G' * W * yi(:));
	xhat = embed(xhat, mask);
	im(pl+4, abs(xhat), 'xhat'), cbar
prompt
end


%
% iterations
%
if ~isvar('xcg'), printm 'qpwls_pcg1'

	xinit = zeros(size(x));
	niter = 100;

	xcg = qpwls_pcg1(xinit(mask), G, W, yi(:), R.C, 'niter', niter);
	xcg = embed(xcg(:,end), mask);

	im(pl+5, abs(xcg), 'matlab qpwls\_pcg1'), cbar

	% examine matlab convergence to xhat
	if isvar('xhat'), printm 'matlab vs xhat'
		diff = xcg - xhat;
		im(pl+6, diff, 'xcg - xhat'), cbar
	end

	max_percent_diff(xcg, xhat)
end
