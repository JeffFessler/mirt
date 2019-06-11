% ir_ct_roi_split1.m
% "frequency split" approach to ROI recon of Lin Fu et al. from Fully 3D 2015.
% Low frequency part of ROI sinogram comes from reprojecting FBP ROI image,
% whereas high frequency part comes from original sinogram.
% ROI reconstructed using PWLS-OS-LALM
% 2015-06-17 Jeff Fessler, University of Michigan

if ~isvar('Af'), printm 'setup geometry'
	f.down = 4;
	igf = image_geom('nx', 512, 'fov', 50, 'down', f.down);
	igf.mask = igf.circ > 0;
	sg = sino_geom('ge1', 'units', 'cm', 'strip_width', 'd', ...
		'down', f.down);

	igr = igf;
	igr.mask = igr.circ(11, 11, -10, 4) > 0;
	im(igf.mask + igr.mask)

	% system objects
	if has_mex_jf
		Af = Gtomo2_dscmex(sg, igf); % full
		Ar = Gtomo2_dscmex(sg, igr); % roi
	else
		Af = Gtomo_nufft_new(sg, igf);
		Ar = Gtomo_nufft_new(sg, igr);
	end
end


if ~isvar('xtrue'), printm 'xtrue, sinogram'
	% read image
	ddir = path_find_dir([filesep 'data']);
	xtrue256 = fld_read([ddir filesep 'ncat,256,slice,140,ct,x100.fld']);
	xtrue256 = single(xtrue256) / 200 * 0.4; % convert to 1/cm units

	if 1 % more realistic sinogram from finer image, avoid "inverse crime"
		ig_big = image_geom('nx', 512, 'fov', igf.fov, 'down', 2);
		if has_mex_jf
			Abig = Gtomo2_dscmex(sg, ig_big);
		else
			Abig = Gtomo_nufft_new(sg, ig_big);
		end
		sino_true = Abig * xtrue256;
	end
	xtrue = downsample2(xtrue256, f.down/2);

	im plc 2 2
	clim = [0 0.4];
	im(1, xtrue, 'x true', clim), cbar
	xlabelf('units: 1 / %s', sg.units)
	im(2, sino_true, 'sino true'), cbar
	im(3, xtrue + igf.mask + igr.mask)

	clear ddir ig_big Abig
drawnow
end


mask2 = conv2(single(igr.mask), ones(9)/9^2, 'same') > 0.999; % for roi rmse
im(3, igf.mask + igr.mask + mask2, 'ROIs')
xl = @(x) xlabelf('RMSE = %.3f / %s', rms(col(x - xtrue)), sg.units);
xr = @(x) xlabelf('RMSE = %.3f / %s', rms(x(mask2) - xtrue(mask2)), sg.units);


if ~isvar('sino'), printm 'noisy sinogram'
	rng(0)
	% transmission data:
	I0 = 1e5;
	yi = poisson(I0 * exp(-sino_true), 0, 'factor', 0.1);
	if any(yi(:) == 0)
		warn('%d of %d values are 0 in sinogram!', ...
			sum(yi(:)==0), length(yi(:)));
	end
	sino = log(I0 ./ max(yi,1)); % noisy fan-beam sinogram
	im(4, sino, 'sino noisy'), cbar
drawnow
%	ir_savefig ir_ct_roi_split1a
end


if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, igf);
	fbp = fbp2(sino, tmp, 'window', 'hanning,0.75');
	im(2, fbp, 'FBP Hanning', clim), cbar
	xl(fbp)
prompt
end


if ~isvar('kappa'), printm 'kappa: try to make resolution approximately uniform'
	wi = yi; % will give 0 weight to any ray where yi=0!
	kappa = sqrt( div0(Af' * wi, Af' * ones(size(wi))) );
	im(3, kappa, 'kappa'), cbar
prompt
end


% use local psf to help select beta
if ~isvar('R'), printm 'R'
	f.l2b = 10; % maybe a bit too big, but ok for now
	f.delta = 0.001;
%	f.pot_arg = {'lange3', f.delta}; % todo: why not as sharp as hyper3?
	f.pot_arg = {'hyper3', f.delta};
	R = Reg1(kappa, 'beta', 2^f.l2b, 'pot_arg', f.pot_arg);
	Rr = Reg1(kappa .* igr.mask, 'beta', 2^f.l2b, 'pot_arg', f.pot_arg);
%	qpwls_psf(A, R, 1, igf.mask, Gdiag(wi), 'loop', 1); % use this to choose beta
end

f.niter = 6;
f.nblock = 41; % 41 subsets

% OS-LALM
if ~isvar('xlalmf'), printm 'iterative reconstruction - lalm full'
	Ab = Gblock(Af, f.nblock);
	xlalmf = ir_pwls_os_lalm(fbp(igf.mask), Ab, sino, R, 'wi', wi, ...
		'isave', 'last', 'niter', f.niter);
	xlalmf = igf.embed(xlalmf);
	im(4, xlalmf(:,:,end), 'LALM full', clim), cbar
	xl(xlalmf(:,:,end))
%	ir_savefig ir_ct_roi_split1b
end


if ~isvar('fbpr'), printm 'fbp 2d fan-beam reconstruction - roi'
	tmp = fbp2(sg, igr);
	fbpr = fbp2(sino, tmp);
	im(1, xtrue .* igr.mask, 'True roi', clim), cbar
	im(2, fbpr, 'FBP roi', clim), cbar
	xr(fbpr)
prompt
end

if ~isvar('Hhi'), printm 'filters'
	nf = 2 * sg.nb;
	u = [-nf/2:nf/2-1]'/nf;
	cut = 0.1;
	Hhi = min((u/cut).^2, 1);
	Hlo = 1 - Hhi;
	plot(u, Hlo, '-o')
	clear u
end

if ~isvar('sino_roi'), printm 'sino_roi'
	sino_f = fft(sino, nf, 1);
	sino_hi = ifft(sino_f .* repmat(ifftshift(Hhi), [1 sg.na]), [], 1);
	sino_hi = sino_hi(1:sg.nb, :);
	im(1, fftshift(sino_f, 1))
	im(2, sino_hi, 'Hi-pass sino')

	sino_Ar = Ar * fbpr;
	sino_Arf = fft(sino_Ar, nf, 1);
	im(1, sino_Ar, 'Reproj FBP')
        sino_lo = ifft(sino_Arf .* repmat(ifftshift(Hlo), [1 sg.na]), [], 1);
	sino_lo = sino_lo(1:sg.nb, :);
	im(3, sino_lo, 'Low-pass sino')

	sino_roi = sino_lo + sino_hi;
	im(4, sino_roi, 'Synth ROI sino')
	clear sino_f sino_Arf
%	ir_savefig ir_ct_roi_split1c
prompt
end


% OS-LALM
if ~isvar('xlalmr'), printm 'iterative reconstruction - lalm roi'
	Ab = Gblock(Ar, f.nblock);
	xlalmr = ir_pwls_os_lalm(fbpr(igr.mask), Ab, sino_roi, Rr, 'wi', wi, ...
		'isave', 'last', 'niter', 2*f.niter);
	xlalmr = igr.embed(xlalmr);
end

if 1 % pics
	im(4, xlalmr(:,:,end), 'PWLS-OS-LALM roi', clim), cbar
	xr(xlalmr(:,:,end))

	im(1, xtrue .* igr.mask, 'True roi', clim), cbar
	im(2, fbpr, 'FBP Ramp roi', clim), cbar
	xr(fbpr)
	im(3, fbp .* igr.mask, 'FBP Hanning roi', clim), cbar
	xr(fbp)
%	ir_savefig ir_ct_roi_split1d
end
