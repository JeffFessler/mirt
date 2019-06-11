% homodyne_recon_demo.m
%
% - Demonstrates use of homodyne_recon.m, which implements the homodyne
% reconstruction method outlined in Doug Noll's 1991 IEEE T-MI paper
% "Homodyne detection in magnetic resonance imaging"
% - Figure shows the decrease in RMSE that results from sampling more.
% - Reconstructed image can be rectangular, with odd or even dimensions.
%
% 2012-06-15, Mai Le
% 2012-09-26 Jeff Fessler tweaked

% set up "true" brain image
%f.dir = [path_find_dir('mri') '/../data/mri/'];
%f.xtrue = [f.dir 'brainweb_t1.jpg'];
%mag_true = double(imread(f.xtrue))'; % true image magnitude
mag_true = 25*ellipse_im(256);
mag_true = mag_true(16:end-16,3:end-2); % make it rectangle to stress test
dims = size(mag_true);

% simulate linear phase over object
[xx yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
ph_max = pi/2;
ph = ph_max*xx/(dims(1)/2) + ph_max*yy/(dims(2)/2);
image = mag_true .* exp(1i*ph); % modulate to make it complex
im_fft = fftshift(fft2(image));

im plc 2 3
im(1, mag_true)

for homodyne_direction = [1 2] % partial k-space can be in direction 1 or 2
	% # of samples in "half" kspace:
	nhalf = ceil((dims(homodyne_direction)+1)/2);
	overlaps = dims(homodyne_direction) - nhalf; % maximum # of extra rows
	overlaps = [0:5:(overlaps-1) overlaps];
	nover = numel(overlaps);

	nrms_im = zeros(nover,1);
	nrms_ph_demod = zeros(nover,1);

	for ii = 1:nover
		overlap = overlaps(ii);
		if (homodyne_direction == 1)
			k1 = 1:min(dims(1), nhalf+overlap);
			partial_kspace = im_fft(k1,:);
		elseif (homodyne_direction == 2)
			k2 = 1:min(dims(2), nhalf+overlap);
			partial_kspace = im_fft(:,k2);
		else
			fail 'bug'
		end

		f.sigma = 800; % complex AWGN std dev
		noisy_partial_kspace = partial_kspace + f.sigma * ...
			(randn(size(partial_kspace)) + 1i*randn(size(partial_kspace)));

		[recon_ph_demod, recon_im, lp_im, full_kspace] = ...
			homodyne_recon(noisy_partial_kspace, dims(1), dims(2), ...
				overlap, 'direction', homodyne_direction);

		nrms_fun = @(x,y) 100 * nrms(x(:), y(:));
		nrms_im(ii) = nrms_fun(abs(recon_im), mag_true);
		nrms_ph_demod(ii) = nrms_fun(abs(recon_ph_demod), mag_true);

		if overlap <= 15
			clim = minmax(mag_true)';
			im(2, abs(recon_im), clim)
			titlef('Conventional, overlap=%d', overlap)
			im(3, abs(recon_ph_demod), clim)
			titlef('Demodulated, overlap=%d', overlap)
			im(4, log(abs(partial_kspace)), clim)
			im(5, abs(recon_im) - mag_true, [-9 9], 'Error')
			im(6, abs(recon_ph_demod) - mag_true, [-9 9], 'Error')
			drawnow
	%		prompt
		end
	end

	f.snr = sqrt(mean(abs(partial_kspace(:)).^2)) / ...
		sqrt(mean(abs(noisy_partial_kspace(:) - partial_kspace(:)).^2));
	f.snr = 20*log10(f.snr);
	pr f.snr

	im subplot 4
	plot(overlaps, nrms_im, 'g-+', overlaps, nrms_ph_demod, 'b-o')
	axis tight
	xtick([0 15 max(overlaps)])
	title('NRMSE (%)')
	legend('without demodulation', 'phase demodulated');
	xlabel('rows acquired beyond 1/2 of kspace');
	ylabel('NRMS error of reconstructed image (%)');
prompt
end
