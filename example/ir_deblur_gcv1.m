% ir_deblur_gcv1.m
% Edge-preserving regularized image deblurring (restoration)
% using (nonlinear!) generalized cross validation GCV (NGCV)
% to select the regularization parameter automatically
% 2013-10-01 Jeff Fessler, University of Michigan

if ~isvar('xtrue')
	ig = image_geom('nx', 100, 'ny', 128, 'dx', 2, 'down', 2);
	xell = [0 0 85 115 0 100;
		0 -60 30 20 0 -80; % mouth
		30 20 25 35 30 20; % eyes
		-30 20 25 35 -30 20;
		35 25 7 7 0 -100; % pupils
		-15 25 7 7 0 -100;
		0 75 60 15 0 -50; % hat
		]; % it is not north park, it is ...
	xtrue = ellipse_im(ig, xell, 'oversample', 3);
	xtrue = single(xtrue);
end

if ~isvar('y')
	psf = ones(3)/9;
	A = Gblur(ig.mask, 'psf', psf);
	ytrue = A * xtrue;

	%% add noise
	rng(0)
	snr = 20; % specify desired SNR of data in dB
	sig = 10^(-snr/20) * norm(ytrue(:)) / sqrt(numel(ytrue));
	y = ytrue + sig * randn(size(ytrue));

	snr_fun = @(x, xtrue) 20*log10(norm(xtrue(:)) / norm(x(:) - xtrue(:)));
	pr 'snr_fun(y, ytrue)'
	pr 'snr_fun(y, xtrue)'

	[nx ny] = size(y);
	nd = nx * ny; % # of data points
%	C = Cdiffs([], 'mask', ig.mask, 'offsets', [1 nx], 'type_diff', 'circshift');

	im plc 2 3
	clim = [-10 130]; % display all images on same gray scale
	im(xtrue, clim)
	im(ytrue, clim)
	im(y, clim)
prompt
end


delta = 5; % small compared to max(y) - min(y)
%dpot = @(t) t ./ sqrt(1 + (t / delta).^2); % derivative of potential
niter = 200; % max # of iter
tol = 1e-5; % stop tolerance for change in x

%% try several values of the regularization parameter

sig = 1; % pretend we don't know the noise level, using sigma=1

reglist = 2 .^ [-3:0.5:2];
nr = numel(reglist);
gcv = nan(nr,1);
rss = nan(nr,1);
snr = nan(nr,1); % solver only
gcv_best = inf;
x_best = [];
ir_best = [];

% emperically it seems better to generate w once outside the loop
% instead of making a new random bernoulli for every beta
w = 2 * (rand(nx,ny) > 0.5) - 1; % bernoulli +/- 1
weps = 0.05 * norm(y(:)) / norm(w(:));

for ir=1:nr
	reg = reglist(ir);
	step = sum(psf(:)) / (1 + reg * 8);

	R = Reg1(ig.mask, 'offsets', [1 nx], ...
		'beta', reg, 'pot_arg', {'hyper3', delta});

	% first run of PCG
	xinit = y;
	x = pwls_pcg1(xinit(ig.mask), A, 1, y(:), R, ...
		'niter', niter, 'stop_diff_tol', tol);
	x = ig.embed(x);

	snr(ir) = snr_fun(x, xtrue);

	% second run of PCG, now with a random perturbation
	xe = pwls_pcg1(xinit(ig.mask), A, 1, col(y + weps * w), R, ...
		'niter', niter, 'stop_diff_tol', tol);
	xe = ig.embed(xe);

	% todo: find reference for thie NGCV formula
	tr_ngcv = (1/nd) * ir_dot_double(w, A * (xe - x) / weps); % NGCV

	rss(ir) = norm(col(A*x) - y(:))^2 / sig^2; % RSS
	gcv(ir) = rss(ir) / (1 - tr_ngcv)^2;

	if gcv(ir) < gcv_best
		gcv_ir_best = ir;
		gcv_best = gcv(ir);
		x_best = x;
	end


	im(x, clim)
	xlabelf('$\log_2(\beta)$ = %4.1f', log2(reg))

%	im subplot 5
%	plot(log2(reglist), rss, '-o')
%	titlef 'RSS', xlabelf '$\log_2(\beta)$'

	im subplot 5
	plot(log2(reglist), snr, '-o', ...
		log2(reglist(gcv_ir_best)), snr(gcv_ir_best), '*')
	titlef 'SNR', xlabelf '$\log_2(\beta)$'

	im subplot 6
	plot(log2(reglist), gcv, '-o', ...
		log2(reglist(gcv_ir_best)), gcv(gcv_ir_best), '*')
	titlef 'NGCV', xlabelf '$\log_2(\beta)$'

	drawnow
end

	im(4, x_best, clim)
	xlabelf('$\log_2(\beta)$ = %4.1f', log2(reglist(gcv_ir_best)))
