% df_example1.m
% parallel-beam image reconstruction:
% comparing direct Fourier (DF) reconstruction vs FBP.
% Copyright 2005-8-10, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'simulate data'
	down = 4;
	% slightly larger than 80cm fov
	sg = sino_geom('par', 'nb', 1500, 'na', 700, 'dr', 0.55, ...
		'strip_width', 'd', 'down', down);
	ig = image_geom('nx', 512, 'dx', [], 'fov', 800, 'down', down);
	ig.dy = ig.dx; % DF needs this for now

	f.erot = 0;
%	ell = []; f.erot = 90;
	ell = []; f.erot = 0;
%	ell = [80 160 50 200 20 1];
%	ell = [0 0 50 200 20 1]; % symmetric, so real spectrum...
	[xtrue ell] = ellipse_im(ig, ell, 'rot', f.erot, 'oversample', 2);

	sino = ellipse_sino(sg, ell, 'oversample', 4);

	im plc 2 2
	clim = [0.9 1.1]; elim = [-1 1]*0.1;
	im(1, ig.x, ig.y, xtrue, 'True', clim), cbar
	im(3, sg.s, sg.ad, sino, 'Sinogram'), cbar
prompt
end

% fbp reconstruction
if ~isvar('xfbp'), printm 'fbp'
	tmp = fbp2(sg, ig);
	cpu etic
	xfbp = fbp2(sino, tmp);
	cpu etoc 'fbp time'

	xrange = @(x) xlabelf('[%g %g]', min(x(:)), max(x(:)));
	im(2, ig.x, ig.y, xfbp, 'FBP', clim), cbar
	im(4, ig.x, ig.y, xfbp-xtrue, 'FBP error', elim), cbar
	xrange(xfbp-xtrue)
%	ir_savefig fig_df1
prompt
end

if 1, printm 'df' % DF recon
	if ir_is_octave
		ilist = {'linear', 'cubic'};
	else
		ilist = {'*linear', '*cubic'};
	end
	olist = {1, 2, 3}; % over-sampling factors
	im plc 3 3
	for ii=1:length(ilist)
	 for oo=1:length(olist)
		ex = {'over', olist{oo}, 'interp', ilist{ii}};
		tmp = fbp2(sg, ig, 'type', 'df,pull', 'extra', ex);
		cpu etic
		xd = fbp2(sino, tmp);
		tim = cpu('etoc');

		if im
			t = sprintf('DF over=%d %s', olist{oo}, ilist{ii});
			im(ii+(oo-1)*3, ig.x, ig.y, xd, t, clim), cbar
			t = sprintf('time %3.1f', tim);
			xlabel(t), drawnow
		end
%		im(6, x, y, xd-xtrue, 'DF error', elim), cbar
%		xrange(xd-xtrue);
	 end
	end
%	ir_savefig fig_df2
return
end


if 1, printm 'df' % explore DF recon
	nu = nx; nv = ny;
%	nu = nx/16; nv = ny/16; % for debug
	nu = 2*nx; nv = 2*ny; % try over-sampling k-space ? (did not help)
	% trick: go to n/2 for simplicity later
	u = [-nu/2:nu/2]/nu/dx;
	v = [0:nv/2]/nv/dy;
	[uu vv] = ndgrid(u, v);
%	npad = nb;
	npad = 2^ceil(log2(2*nb-1))
	rho = [-npad/2:npad/2-1]'/npad/dr;
	phi = deg2rad(f.orbit_start + [0:na-1]/na * f.orbit);
	qrho = sqrt(uu.^2 + vv.^2);
	qphi = atan2(vv, uu);
%	im clf, plot(qrho, qphi, '.', max(rho), pi/2, 'x'), return % see samples
%	minmax(qphi)

	ws = (nb-1)/2 + f.offset_s;
%	trick: because ws is non-integer, we must fftshift the phase as follows:
	phase1 = dr * fftshift(exp(2i*pi*[-npad/2:npad/2-1]'/npad * ws));

	% phase corresponding to image-domain half-pixel shift 
	f.offset_x = 0;
	f.offset_y = 0;
	wx = (nx-1)/2 + f.offset_x;
	wy = (ny-1)/2 + f.offset_y;
	phase2 = outer_sum([-nu/2:nu/2]/nu * (wx-nx/2), [0:nv/2]/nv * (wy-ny/2));
	phase2 = exp(-2i*pi*phase2);

	cpu etic
	f_sino = fft(sino, npad, 1);
	f_sino = repmat(phase1, [1 na]) .* f_sino;
%	plot(angle(f_sino(:,1))), return
%	f_sino = reale(f_sino, 1e-12); % for testing only
	f_sino = fftshift(f_sino, 1);

	% trick: add another row at pi by flipping, to help interpolator
	if 1
		if f.orbit ~= 180 || f.orbit_start ~= 0, error 'need 0,180', end
		f_sino = [f_sino, [conj(f_sino(1,1)); flipud(f_sino(2:end,1))]];
		phi = [phi pi];
	end
%	im clf, im(rho, phi, log(abs(f_sino))), cbar
%	im(3, rho, phi, abs(f_sino), '|Sino|'), cbar
%	im(3, rho, phi, real(f_sino), 'real(Sino)'), cbar
%	im(6, rho, phi, imag(f_sino)), cbar, return

	f.interp = '*nearest';
	f.interp = '*spline';
	f.interp = '*linear';
	f.interp = '*cubic';
	% [-nu/2,nu/2] x [0,nv/2]
	f_xd = interp2(phi, rho, f_sino, qphi, qrho, f.interp);
	if any(isnan(f_xd)), error 'nan', end

	f_xd = phase2 .* f_xd;

	if 1 % edges of spectrum must be real for iDFT to be real
		f_xd(1+[0:2]*nu/2, nv/2+1) = real(f_xd(1+[0:2]*nu/2, nv/2+1));
		f_xd(1+[0:2]*nu/2, 1) = real(f_xd(1+[0:2]*nu/2, 1));
		f_xd(1:nu/2,nv/2+1) = conj(flipud(f_xd(nu/2+2:end,nv/2+1)));
	end

	% form entire spectrum using conjugate symmetry
	t1 = fliplr(flipud(f_xd(2:end-1, 2:end))); % [-nu/2+1,nu/2-1] x [1,nv/2]
	t1 = [fliplr(f_xd(1, 2:end)); t1]; % [-nu/2] x [1,nv/2] added
	t2 = f_xd(1:end-1, 1:end-1); % [-nu/2,nu/2-1] x [0,nv/2-1]
%	t2 = f_xd(1:end-1, 2:end-1); % [-nu/2,nu/2-1] x [1,nv/2-1]
%	t0 = f_xd(nu/2+1:end, 1); % [0,nu/2] x [0]
%	t0 = [conj(flipud(t0(2:nu/2+1))); t0(1:nu/2)]; % [-nu/2,-1] x [0] added
	f_xd = [conj(t1), t2];

%	im clf, dft_sym_check(fftshift(f_xd)), return
%	im clf, im(f_xd), return

	if 0 % examine errors
		f_true = fftshift(fft2(fftshift(xtrue))) * dx * dy;
%		f_xd = fftshift(fft2(fftshift(xfbp))) * dx * dy;
		im(131, log(abs(f_true)), 'f true'), cbar h
		im(132, log(abs(f_xd)), 'f xd'), cbar h
		im(133, log(abs(f_xd - f_true)), 'f err'), cbar h
%		[iu iv] = imax(abs(f_xd - f_true), 2);
		iu = 1:nu; iv = nv/2+0;
		if im
			clf, plot(iu, abs(f_true(iu,iv)), '-o', ...
			iu, abs(f_xd(iu,iv)), '.')
			axis([230 280 0 3e4])
		end
	return
	end
	if 0
		t1 = f_xd(2:end,2:end);
		t2 = conj(flipud(fliplr(t1)));
		im clf, im([t1 t2 t1-t2]')
	return
	end

	xd = fftshift(ifft2(fftshift(f_xd))) / dx / abs(dy);
	xd = reale(xd, 'warn');
	if (nu == 2*nx)
		xd = xd([0:nx-1]+1+nx/2, [0:ny-1]+1+ny/2);
	end
	cpu etoc 'df time'

	max_percent_diff(xd, xd2)

	im(3, x, y, xd, 'DF', clim), cbar
	im(6, x, y, xd-xtrue, 'DF error', elim), cbar
	xrange(xd-xtrue);

	nrms(xtrue, xfbp)
	nrms(xtrue, xd)
return
end

% df reconstruction via nufft
if 0 || ~isvar('xdf'), printm 'dfr'
	mask = true(nx,ny);
	f.kb_m = [2 2];
	f.kb_alf = 14.0 * [1 1]; % fix: need to tune!
	f.interp = {'table', 2^11, 'kaiser', f.kb_alf, f.kb_m};
%	f.interp = {'table', 2^11, 'kaiser'};
%	f.interp = {'table', 2^11, 'minmax:kb'}; % for forward projection
	Gn = Gtomo_nufft(mask, [nb na], ...
		'J', 2, ...
		'pixel_size', dx, 'ray_spacing', dr, ...
		'orbit', f.orbit, 'orbit_start', f.orbit_start, ...
		'interp', f.interp);
	if im
		im subplot 1
		plot([real(Gn.arg.st.h{1}) imag(Gn.arg.st.h{1})])
		im(2, Gn.arg.st.sn), cbar
	end

	if 0 % check that forward projection matches
		s2 = Gn * xtrue;
		im([s2 - sino])
	end

	% density compensation
	om = Gn.arg.st.om;
%	wt = ir_mri_density_comp(om, 'voronoi');
%	wt = reshape(wt, [], na);
	om = reshape(om, [], na, 2);
%	im clf, plot(om(:,:,1), om(:,:,2), '.')
	wt = sqrt(om(:,:,1).^2 + om(:,:,2).^2);
%	im(wt)
%	im(wt)
	Gn.arg.tomo_filter = wt .* Gn.arg.tomo_filter;

	cpu etic
profile on
	xdf = Gn' * sino; % no! need to apply density compensation
profile report
	cpu etoc 'df time'
%	xdf = mean(xtrue(:)) / mean(xdf(:)) * xdf;
	im(3, xdf, 'df', clim), cbar
	im(3, xdf, 'df'), cbar
	im(6, xdf-xtrue, 'df error'), cbar
prompt
end

return

if ~has_aspire, return, end % check consistency with aspire

if ~isvar('im_asp')
	dir = test_dir;
	f.sino = [dir 'sino.fld'];
	f.image = [dir 'image.fld'];
	f.dsc = [dir 't.dsc'];
	fld_write(f.sino, sino)

	t = geom.arg_dsc;
	t(1,8) = '2'; % trick: replace system 9 with system 2

	char_array_write(f.dsc, t)
	f.win = 'boxcar,1,0,1';
	com = sprintf('echo y | i -chat 0 fbp dsc %s %s %s %s', ...
		f.image, f.sino, f.dsc, f.win);
	cpu etic
	disp(os_run(com))
	cpu etoc 'aspire time'

	im_asp = fld_read(f.image);
	im(1, im_asp, 'aspire', clim), cbar
	max_percent_diff(xfbp, im_asp)
prompt
end

if im
	iy = ny/2; ix = 1:nx;
	subplot(133)
	plot(ix, xtrue(ix,iy), '-', ix, recon(ix,iy), '--')
	axis([1 nx -0.5 10.5]), legend('true', 'recon', 4)
end
