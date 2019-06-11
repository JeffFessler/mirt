% mri_r2_fit_ist1.m
%| Fit R2=1/T2 and "pd" image using ODWT regularization
%| and iterative soft thresholding (IST), with alternating updates.
%| Use interactive sliders to select soft threshold value.
%| Copyright 2012-05-20, Jeff Fessler, University of Michigan

%| in
%|	mask
%|	yi	[(N) nt]	data
%|	telist	[nt]		echo times
%|	pdconv	[(N)]		initial guess
%|	r2conv	[(N)]		""

	im plc 2 3
	im(1, yi)
	im(4, r2true)


if ~isvar('U'), printm 'setup ODWT object'
	f.level = 1; % v1
	f.level = 2; % v2
	f.wname = 'haar'; % v1,2
%	f.wname = 'sym2'; % v3
	U = Godwt1(mask, 'level', f.level, 'wname', f.wname)';
	detail = rem(U.codes,10) ~= 0; % which are the "detail" coefficients
%	im(detail)
end

if 0 % show wavelet coefficients (cheat: helps determine scale)
	im plc 1 2
	im(1, U' * r2true), cbar
	im(2, U' * pdtrue), cbar
return
end


	climpd = [0 1];

	f.delta = 0.01;
	pot = potential_fun('fair-l1', f.delta);
%	shrink = @(t, a) pot.shrink(t, a); % v1
	shrink = @(t, a) pot.shrink(t, a .* detail); % v2 shrink only details

	% todo: assess rms error of r2 only within some mask!
	good = pdtrue > 0.1;
	rms_pd = @(pd) xlabelf('RMS = %.2f', rms(pd(good) - pdtrue(good)));
	rms_r2 = @(r2) xlabelf('RMS = %.2f', rms(r2(good) - r2true(good)));

if 0 % explore
	thresh = 2^-4;
	im plc 2 3
	tmp = U * shrink(U' * pdtrue, thresh);
	im(1, pdtrue, climpd); im(2, tmp, climpd); im(3, tmp-pdtrue)
	rms_pd(tmp)
	tmp = U * shrink(U' * r2true, thresh);
	im(4, r2true, clim); im(5, tmp, clim); im(6, tmp-r2true)
	rms_r2(tmp)
return
end


if 1 || ~isvar('xi'), printm 'IST reconstruction'

	pdn = pdconv;
	r2n = r2conv;

	im(2, abs(pdn), climpd), rms_pd(pdn)
	im(3, r2n, clim), rms_r2(r2n)

	if ~isvar('sl_save1')
		sl_save1 = 0.2338; % v1
		sl_save2 = 0.5216;
		sl_save1 = 0.7696; % v2
		sl_save2 = 0.3849;
	end

	h1 = jf_add_slider('pos', [0.1 0.0 0.4 0.03], ...
		'callback', [], 'value', sl_save1);
	h2 = jf_add_slider('pos', [0.5 0.0 0.4 0.03], ...
		'callback', [], 'value', sl_save2);

	hb = uicontrol('style', 'togglebutton', 'string', 'stop',...
		'units', 'normalized', 'position', [0.01 0.0 0.08 0.03]);
	drawnow

%	tfun1 = @(x) 2^-9 + 2^-6 * x;
	tfun1 = @(x) 2^-11 + 2^-7 * x; % for detail = 1
	tfun2 = @(x) 2^-8 + 2^-3 * x;


	curv1 = numel(telist); % because exp(-*) < 1
	curv2 = sum(telist.^2); % because exp(-*) < 1 and |pd| <= 1 approx

	f.show = 5;
	iter = 0;
	while (iter < 1000)
		if (get(hb, 'value')) % stop button
			break
		end

		sl1 = get(h1, 'value');
		sl2 = get(h2, 'value');
		if sl1 ~= sl_save1 ...
		|| sl2 ~= sl_save2
			sl_save1 = sl1;
			sl_save2 = sl2;
			iter = 0;
		end
		thresh1 = tfun1( sl1 );
		thresh2 = tfun2( sl2 );

		% todo: try FISTA instead of IST

		% IST update for pd
		tmp = struct('yi', yi, 'r2', r2n, 'te', telist); % data
		[dummy tmp] = mri_r2_fit_costgrad_pd(pdn, tmp); % grad
		tmp = pdn - 1/curv1 * tmp;
		pdn = U * shrink(U' * tmp, thresh1);

		% IST update for r2
		tmp = struct('yi', yi, 'pd', pdn, 'te', telist); % data
		[dummy tmp] = mri_r2_fit_costgrad_r2(r2n, tmp); % grad
		tmp = r2n - 1/curv2 * tmp;
		r2n = U * shrink(U' * tmp, thresh2);

		if ~rem(iter, f.show)
			im(5, abs(pdn), climpd), rms_pd(pdn)
			ylabelf('iter = %d', iter)
			titlef('thresh = %g', thresh1)

			im(6, r2n, clim), rms_r2(r2n)
%			ylabelf('iter = %d', iter)
			titlef('thresh = %g', thresh2)
			drawnow
		end

		if 0 && ~rem(iter, f.show)
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
%	xi = embed(xi, mask);
end

if 1
	im plc 1 3
	im(1, r2true, clim, 'True R2 map')
	im(2, r2conv, clim, 'Conventional R2 fit'), rms_r2(r2conv)
	im(3, r2n, clim, 'Regularized R2 estimate'), rms_r2(r2n)
	tmp = strrep(num2str(telist*1e3), '  ', '_');
	tmp = sprintf('fit_t2_v1_r2_snr=%d_te=%s', f.snr, tmp)
%	ir_savefig(tmp)
end

if 0
	im clf
	im('row', 1, yi)
	titlef('Noisy spin-echo images for TE = %s ms', num2str(telist*1e3))
	tmp = strrep(num2str(telist*1e3), '  ', '_');
	tmp = sprintf('fit_t2_v1_yi_snr=%d_te=%s', f.snr, tmp)
%	ir_savefig('cw', tmp)
end
