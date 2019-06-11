% ir_inpaint_admm1_test
% illustrate using "ADMM" for image inpainting (aka Split Bregman algorithm)
% uses l1 potential function by default, i.e., anisotropic total variation (TV)
% The cost function here is
% $cost(x) = 1/2 \| y - diag(samp) x \|_2^2 + \beta * pot(C x)$

if ~isvar('yi'), printm 'generate noisy image'
	tmp = 'neuhoff.jpg';
	if exist(tmp, 'file') && 1
		xt = single(imread(tmp))';
	else
		xt = 250/9 * ellipse_im(200);
	end
	clim = [0 255];
	im(xt, clim)
	titlef 'Original image'
%	ir_savefig fig_xt

	snr = @(x) -20*log(norm(x(:)-xt(:)) / norm(xt));
	xlsnr = @(x) xlabelf('SNR = %4.1f dB', snr(x));

	siz = size(xt);
	samp = false(siz);
	samp(2:4:end,:) = true;
	samp(:,2:4:end) = true;
	im(samp)
	titlef 'Cutset sampling pattern'
%	ir_savefig fig_samp1

	S = Gdiag(samp);
	Ss = sparse(S);
	Ss = Ss(samp(:),:);
%	yt = samp .* xt;
	yt = S * xt;
	im(yt, clim)

	if 1
		yi = yt; % noiseless
		im(yi, clim)
		titlef 'Cutset sampled image'
		xlabelf '("noiseless")'
%		ir_savefig fig_yt
	else
		rng(7)
		sig = 0;
		yi = yt + sig * randn(size(yt));

		im plc 2 2
		im(1, xt, clim, 'True')
		im(2, yi, clim, 'Noisy')
		xlabelf('SNR = %4.1f dB', snr(yi))
	end

end


if ~isvar('x0'), printm 'initialize'
	[xx yy] = ndgrid(1:siz(1), 1:siz(2));
	x0 = griddata(xx(samp), yy(samp), yi(samp), xx, yy);
	x0(isnan(x0)) = mean(x0(~isnan(x0)));
	clear xx yy
	im(x0, clim)
	titlef 'griddata interpolated'
	xlsnr(x0)
%	ir_savefig fig_x0g
end


if ~isvar('C'), printm 'C'
	nx = siz(1);
	ny = siz(2);
	nn = prod(siz);
%	nx = 4; ny = 5;
	Ch = speye(nn) - circshift(speye(nn), 1); % lazy!
	Cv = speye(nn) - circshift(speye(nn), nx);
	C = [Ch; Cv];
%	im(C')
	if 0
		tmp = C * double(xt(:));
		tmp = reshape(tmp, [siz 2]);
		im(tmp)
	end
end


if ~isvar('x1'), printm 'initialize gmrf'
	reg1 = 2^-3;
	H = Ss'*Ss + reg1 * C'*C;
	x1 = H \ (Ss' * double(yi(samp)));
	x1 = reshape(x1, siz);
	im(x1, clim)
	titlef 'Simple Gaussian MRF reconstruction'
	xlsnr(x1)
%	ir_savefig fig_x1m
end


if 0 || ~isvar('xs'), printm 'run admm'
	[xs snrs] = ir_inpaint_admm1(yi, samp, 'beta', 2^-4, ...
		'x0', x0, ...
		'stop_diff_tol', 1e-4, 'chat', 0, ...
		'isave', 'all', ...
		'userfun', @(x, iter) snr(x), ...
		'niter', 2^7, 'shrink', [], 'rho', 2^-2);
	xh = xs(:,:,end);
	im(xh, clim)
	titlef 'Sparsity-based inpainting with ADMM'
	xlsnr(xh)
%	ir_savefig fig_xh
end


if 1 % figures
	im plc 2 3
	im(1, xt, clim, 'True x')
	im(2, samp, 'Cutset sampling')
	im(3, yi, clim, 'Cutset sampled')
	im(4, x0, clim, 'griddata')
	xlsnr(x0)
	im(5, x1, clim, 'Gaussian MRF')
	xlsnr(x1)
	im(6, xh, clim, 'ADMM')
	xlsnr(xh)

%	ir_savefig fig_all_sl200
return

	im subplot 4
	iters = 0:numel(snrs)-1;
	plot(iters, snrs, '-o')
	axis([0 20 20 55])
	xlabel 'ADMM iteration'
	ylabel 'SNR (dB)'
	ytick([25 53])
%	ir_savefig cw ir_denoise_admm1_test
end


if 0 % movie of iterates
	tmp = xs(:,:,[1:20 end]);
	niter = size(tmp,3);
	bar = size(tmp,2) / niter;
	for ii=1:niter
                tmp(1:3, floor([1:bar]+(ii-1)*bar), ii) = 8;
        end
	clf, im(tmp)
%	movie2(tmp, 'clim', clim, 'fps', 10)
end
