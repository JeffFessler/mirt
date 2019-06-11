% demo_resolution_noise_tradeoff
% Show tradeoff between spatial resolution and noise using
% 2D fan-beam CT example, for original and modified regularization
% Copyright 2010-03-27, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'simulate data'

	ig = image_geom('nx', 2^5, 'fov', 48);
%	ig.mask = ig.circ > 0; % reconstruction support: a disk

	% fan-beam sinogram geometry approximating GE Lightspeed system
	sg = sino_geom('ge1', 'units', 'cm', 'na', 256, 'down', 8);
	sg.strip_width = sg.ds;

	% shepp-logan image and parameters
	ell = [0 0 20 15 0 0.2;
		-9 0 5 5 0 -0.1;
		9 0 5 5 0 0.1];
	xtrue = ellipse_im(ig, ell, 'oversample', 2);

	% analytical sinogram with ray subsampling to model detector size
	sino = ellipse_sino(sg, ell, 'oversample', 4);

	im plc 2 3
	clim = [0 0.22];
	im(1, ig.x, ig.y, xtrue, 'True', clim), cbar
	im(2, ig.x, ig.y, ig.mask, 'Mask'), cbar
	im(4, sg.s, sg.ad, sino, 'Ideal sinogram'), cbar
prompt
end


if ~isvar('A'), printm 'A' % system matrix
	A = Gtomo2_wtmex(sg, ig, 'nthread', jf('ncore'));
	wi = exp(-sino); % data weighting for monoenergetic CT
	W = Gdiag(wi);
end

% standard regularizer
if ~isvar('R1'), printm 'R1'
	R1 = Reg1(ig.mask, 'beta', 1, 'offsets', [1 ig.nx]);
	C1 = R1.C1;
	tmp = full(C1);
	hR1 = tmp' * tmp;
end

if ~isvar('F'), printm 'F'
	F = full(A);
	F = F' * W * F;
end

if 0
	reg = 1;
	psf = (F + reg * hR0) \ (F * ej);
	psf = ig.embed(psf);

	tmp = (F + reg * hR0) \ ej;
	cov = tmp' * F * tmp;
	im(psf)
end

% "certainty factors" from Fessler 1996 T-IP spatial resolution paper
if ~isvar('kappa'), printm 'kappa'
	kappa = sqrt( div0(A' * wi, A' * ones(size(wi))) );
	im(3, kappa, 'kappa'), cbar
end

% regularizer with kappa modification
if ~isvar('R2'), printm 'R2'
	R2 = Reg1(kappa, 'beta', 1, 'offsets', [1 ig.nx]);
	C2 = R2.C;

	tmp = full(C2);
	hR2 = tmp' * tmp;
end

prompt

if im
	xpos = ell([2 3 1],1);
	ej = zeros(ig.np, length(xpos));
	for ii=1:length(xpos);
		ej(:,ii) = col(ig.unitv(ig.nx/2+xpos(ii),ig.ny/2+1));
	end
	Fe = F * ej;
	rfun = @(x) 2^(-4 + 5.1 * x);
	im clf
	plot(0, 0, 'w.');
	hp = gca;
	xlabel 'Blur: FWHM(PSF)'
	ylabel 'Noise: std dev'
	ax = [1 3 0 1.4];
	axis(ax)
	hr = uicontrol('style', 'slider', 'string', 'value', ...
		'units', 'normalized', 'position', [0.1 0.0 0.8 0.03], ...
		'value', 0);
	hb = uicontrol('style', 'togglebutton', 'string', 'stop',...
		'units', 'normalized', 'position', [0.01 0.0 0.08 0.03]);
%	ha = axes('pos', [0.7 0.6 0.2 0.1]);
	drawnow
end

	reg_save = nan;
	h1 = [];
	h2 = [];

while (1)
	if get(hb, 'value') % stop button
		break
	end

	reg1 = rfun(get(hr, 'value'));
	if reg_save ~= reg1
		reg_save = reg1;
	end

	reg2 = 400 * reg1;
	psf1 = (F + reg1 * hR1) \ Fe;
	psf1 = ig.embed(psf1);
	psf2 = (F + reg2 * hR2) \ Fe;
	psf2 = ig.embed(psf2);
	for ii=1:length(xpos)
		fw1(ii) = fwhm2(psf1(:,:,ii));
		fw2(ii) = fwhm2(psf2(:,:,ii));
	end

	tmp = (F + reg1 * hR1) \ ej;
	cov1 = tmp' * F * tmp;
	std1 = sqrt(diag(cov1));
	tmp = (F + reg2 * hR2) \ ej;
	cov2 = tmp' * F * tmp;
	std2 = sqrt(diag(cov2));

	color1 = {'ro', 'go', 'bo'};
	color2 = {'r+', 'g+', 'b+'};
	hold on
	for ii=1:length(xpos)
		plot(hp, fw1(ii), std1(ii), color1{ii})
		plot(hp, fw2(ii), std2(ii), color2{ii})
	end
	if ~isempty(h1), delete(h1), end
	h1 = plot(hp, fw1, std1, 'c:');
	if ~isempty(h2), delete(h2), end
	h2 = plot(hp, fw2, std2, 'y:');
	hold off
%	plot(ha, 0, 0, '.')
%	im(psf2, ' ')
	drawnow
end

if 0 % check nominal resolution
	qpwls_psf(A, R, 1, ig.mask, W);
end
