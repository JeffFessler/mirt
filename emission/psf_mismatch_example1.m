% psf_mismatch_example1.m
%
% 1D example showing the effect of PSF mismatch on ML-EM algorithm
%
% Copyright 2001-8-24, Jeff Fessler, University of Michigan


%
% generate data
%
if ~isvar('yi'), printm 'data'
	n.x = 128;	n.y = 1;
	xtrue = zeros(n.x, n.y) + 1;
	xtrue(n.x/2+1+[-5:5]) = 100;

	normcdf = @(x,mu,sig) 0.5 * erfc(-(x-mu)./sig / sqrt(2));
	ghist = @(x, fwhm) normcdf(x+1/2, 0, fwhm/sqrt(log(256))) ...
			- normcdf(x-1/2, 0, fwhm/sqrt(log(256)));

	nk = 21;
	xk = [-(nk-1)/2 : (nk-1)/2]';

	f.fwhm0 = 5;
	f.fwhm1 = 7;
	f.fwhm2 = 2;

	kern0 = ghist(xk, f.fwhm0);
	printf('kern0 discard %g', 1-sum(kern0))
	kern0 = kern0 / sum(kern0);
	G0 = filtmat('1d', kern0, n.x);

	kern1 = ghist(xk, f.fwhm1);
	printf('kern1 discard %g', 1-sum(kern1))
	kern1 = kern1 / sum(kern1);
	G1 = filtmat('1d', kern1, n.x);

	kern2 = ghist(xk, f.fwhm2);
	printf('kern2 discard %g', 1-sum(kern1))
	kern2 = kern2 / sum(kern2);
	G2 = filtmat('1d', kern2, n.x);

	ri = 5;
	yi = G0 * xtrue + ri;
	if im
		plot(1:n.x, yi, '-o'), xlabel i, ylabel 'y_i', title Data
	end
prompt
end

% uniform initial image
xinit = ones(size(xtrue));

%
% ML-EM iterations
%
if ~isvar('x2'), printm 'ML-EM for various psf models'
	f.niter = 81;
	x0 = eml_em(xinit, G0, yi, 1, ri, 'niter', f.niter, 'isave', 'all');
	x1 = eml_em(xinit, G1, yi, 1, ri, 'niter', f.niter, 'isave', 'all');
	x2 = eml_em(xinit, G2, yi, 1, ri, 'niter', f.niter, 'isave', 'all');
	im plc 2 2
	im(1, x0, 'x0')
	im(2, x1, 'x1')
	im(3, x2, 'x2')
prompt
end

	fw_true = fwhm1(xtrue, 'imid', n.x/2+1);
	fw0 = fwhm1(x0(:,2:end), 'imid', n.x/2+1); % gives warning due to hole
	fw1 = fwhm1(x1(:,2:end), 'imid', n.x/2+1);
	fw2 = fwhm1(x2(:,2:end), 'imid', n.x/2+1);

if 1 && im
	ii = 1:f.niter;
	clf, plot(ii(1:4:end), fw0(1:4:end), 'yo', ...
		ii(1:4:end), fw1(1:4:end), 'cx', ...
		ii(1:4:end), fw2(1:4:end), 'g+', ...
		[min(ii) max(ii)], fw_true * [1 1], 'r--', ...
		ii, fw0, 'y-', ...
		ii, fw1, 'c-', ...
		ii, fw2, 'g-')
	legend(	sprintf('True system PSF, FWHM=%g', f.fwhm0), ...
		sprintf('Model PSF too big, FWHM=%g', f.fwhm1), ...
		sprintf('Model PSF too small, FWHM=%g', f.fwhm2), ...
		'ideal FWHM')
	xlabel 'Iteration of ML-EM algorithm'
	ylabel 'FWHM of reconstructed signal'
end
