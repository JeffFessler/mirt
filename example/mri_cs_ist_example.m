%% mri_cs_ist_example.m
%|
%| Compressed-sensing for MRI via iterative soft thresholding (IST).
%| Interactive demo with slider to select soft threshold value.
%|
%| Copyright 2010-02-23, Jeff Fessler, University of Michigan

%% system matrix
if ~isvar('A'), printm 'setup Gdft object'
	Nd = [1 3/2] * 2^7;
	rng(0)
	samp = rand(Nd) > 0.3;
	mask = true(Nd);
	mask(1) = false; % todo: stress
	A = Gdft('mask', mask, 'samp', samp);

	im plc 2 3
	im(samp), cbar
	titlef('%.1f\%% samples', mean(samp(:))*100)
end


%% Orthogonal discrete wavelet transform ODWT
if ~isvar('U'), printm 'setup ODWT object'
	U = Godwt1(mask, 'level', 3)';
	ic = 3500;
	ic = min(ic,prod(Nd));
	im(embed(U(:,ic), mask)), cbar
	titlef('ODWT column')
end


%% true object
clim = [0 9];
if ~isvar('xt'), printm 'setup object'
	xt = ellipse_im(Nd, [], 'oversample', 2);
	im(abs(xt), '|\x| true', clim), cbar

	% show wavelet coefficients (cheat: helps determine scale)
	im(reshape(U' * xt, Nd), [-15 15]), cbar
	titlef('$\U'' \x$')
end


%% parameters
if 1
	pot = potential_fun('l1');
	f.delta = 0.1;
	pot = potential_fun('fair-l1', f.delta);
%	t = linspace(-2, 2, 101);
%	plot(t, abs(t), '-', t, pot.potk(t), '--')

	shrink = @(t, a) pot.shrink(t, a);
%	shrink = @(t, a) (t - a * sign(t)) .* (abs(t) > a);
	thresh = 2^-4;
end


%% explore shrinkage
if 0
	tmp = U * shrink(U' * xt, thresh);
	im plc 1 3
	im(xt, clim); im(tmp, clim); im(tmp-xt)
return
end


%% noisy data
if ~isvar('yi'), printm 'setup data yi'
	yb = A * xt(mask);
	sig = 30;
	yi = yb + sig * randn(size(yb));
	snr = 20 * log( norm(yb(:)) / norm(yi(:) - yb(:)) );
	im(log(max(abs(fftshift(embed(yi,samp))),1)), [1 10]), cbar
	titlef('k-space data')
	xlabelf('SNR = %0.1f dB', snr)
end


%% iFFT recon
if ~isvar('xf'), printm 'basic inverse FFT reconstruction'
	xf = A' * yi / prod(Nd);
	xf = embed(xf, mask);
	im(abs(xf), '|\x| "zero-pad"', clim), cbar
prompt
end

if ~im, return, end

%% Iterative soft thresholding (IST) recon (todo: replace with FISTA!)
if 1 || ~isvar('xi'), printm 'IST reconstruction'

	im plc 2 1, im(1, abs(xf), clim)
	hu = jf_add_slider('pos', [0.1 0.0 0.8 0.03], 'callback', []);

	hb = uicontrol('style', 'togglebutton', 'string', 'stop',...
		'units', 'normalized', 'position', [0.01 0.0 0.08 0.03]);
	drawnow
prompt

	tfun = @(x) 2^-8 + 2^-1 * x;
	thresh_save = nan;

	% todo: somehow threshold only the detail coefficients
	curv = prod(Nd); % because standard DFT

%	reg = thresh * curv;
	xi = xf(mask);
	f.show = 8;
	iter = 0;
	niter = 100;
	while (iter <= niter)
		if (get(hb, 'value')) % stop button
			break
		end
		thresh = tfun( get( hu, 'value') );
		if thresh_save ~= thresh
			thresh_save = thresh;
			iter = 1;
		end
		tmp = xi + 1/curv * A' * (yi - A * xi);
		xi = U * shrink(U' * tmp, thresh);
		if ~rem(iter, 4)
			% thresholded truth for comparison
			tmp = U * shrink(U' * xt, thresh);
			tmp = [tmp; embed(abs(xi), mask)];
			im(1, tmp, clim, ' '), cbar
			fun = @(p,s,v) ...
			text(p*Nd(1)/2, -4, sprintf(s, v), 'horiz', 'cent');
			fun(1, 'True: threshold = %g', thresh)
			fun(3, 'IST %4d', iter)

			wl = [reshape(U' * xt(mask), Nd); ...
				reshape(U' * xi(:), Nd)];
			im(2, abs(wl), [0 4])
			drawnow
		end
		iter = iter + 1;
	end
	xi = embed(xi, mask);
	printm('done')
end


if 0
	im plc 1 3
	im(xt)
	im(xf)
	im(xi)
end
