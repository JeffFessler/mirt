% mri_example_epi_b0.m
% Example of B0 field-corrected MR reconstruction using EPI trajectory.
% WORK IN PROGRESS!
% Copyright 2010-03-10, Jeff Fessler, University of Michigan

if ~isvar('xtrue'), printm 'xtrue'
	f.dir_xtrue = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir_xtrue 'brainweb_t1.jpg'];
	xtrue = single(imread(f.xtrue))' / 230 * 8;
	xtrue = interp2(-129:128, [-129:128]', xtrue, ...
		[-64:63]*1.8, [-64:63]'*1.8);

	[nx ny] = size(xtrue);
	ig = image_geom('nx', nx, 'ny', nx, 'fov', 21); % cm
	ig.mask = ellipse_im(ig, [0 0 10 10 0 1]) > 0;
	pr sum(xtrue(ig.mask) ~= 0)
	xtrue = xtrue .* ig.mask;

	f.clim = [0 9];
	if 1 && im
		im clf, im(xtrue + 5*ig.mask_outline, f.clim), cbar
		title('Image and support outline')
	prompt, clear t
	end
end


if ~isvar('zmap'), printm 'zmap'
	f.dir_fmap = [path_find_dir('mri') '/phase-data/'];
	f.fmap = [f.dir_fmap 'fieldmap128.fld'];
	fmap = fld_read(f.fmap); % [-40 128] Hz
%	im(fmap + xtrue*10)
	zmap = 0 + (2i*pi) * fmap;

	if 1 && im % picture of field map
		im clf, im(fmap + max(fmap(:)) * ig.mask_outline)
		titlef 'Fieldmap and mask outline'
		cbar 'Hz'
	end
prompt
end


% kspace trajectory
if ~isvar('kspace'), printm 'kspace (slow due to voronoi)'
	N = [nx ny];
	f.traj = 'epi-under'; f.dens = {'voronoi'};
	rng(0)
%	samp = rand(ig.ny,1) < 0.8;
	samp = true(ig.ny,1); % fully sampled
	printm('%% samples used: %g', sum(samp) / length(samp) * 100)
	[kspace, omega, wi_traj] = mri_trajectory(f.traj, {'samp', samp}, N, ig.fov);
	if 1 && im
		clf, h = plot(omega(:,1), omega(:,2), '.-');
		set(h, 'markersize', 5), clear h
		axis(pi*[-1 1 -1 1] * 1.05), axis_pipi, axis square
	prompt
	end
end

% sample times, starting at 0
if ~isvar('ti')
	f.dt = 5e-6;
%	f.dt = 0; % ideal
	f.nshot = 3 % todo: improve realism
	ti = single(([1:length(kspace)]-1)*f.dt) / f.nshot; % cheap multi-shot
	f.daq = max(ti) + f.dt;
	printm('nshot=%d, readout time: %g ms', f.nshot, 1000 * f.daq)
end


% "exact" system for generating the data
if ~isvar('Ge_zmap'), printm 'Ge_zmap'
	Ge_ft = Gmri(kspace, ig.mask, 'exact', 1, 'n_shift', N/2, ...
		'fov', ig.fov, 'basis', {'sinc'});
	% trick! adjust wi's to undo the basis effect for CP
	wi_basis = wi_traj ./ Ge_ft.arg.basis.transform;

	Ge_zmap = feval(Ge_ft.arg.new_zmap, Ge_ft, ti, zmap, {});
	clear yi
end


% base NUFFT (no field map)
% todo: why do i need 6x6 hood here?  samples lie on cartesian grid!?
if ~isvar('Gn'), printm 'Gn'
	f.nufft = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
	Gn = Gmri(kspace, ig.mask, 'fov', ig.fov, 'nufft', f.nufft, ...
 		'basis', {Ge_ft.arg.basis.type});

	if 0 % check Gn vs exact
		% todo: does not match for J = 2 or J = 1 !?
		tmp = ig.unitv('c', [0 2]);
		plot([real(Gn * tmp) real(Ge_ft * tmp)])
	end
end

% system
if ~isvar('Gm'), printm 'Gmri'
	L = 9;
	Gm = feval(Gn.arg.new_zmap, Gn, ti, zmap, L);
end


% data
if ~isvar('yi'), printm 'yi'
	yt = Ge_zmap * xtrue(ig.mask);
	% add noise
	rng(0)
	yn = randn(size(yt)) + 1i * randn(size(yt));
	f.snr_db = 50;
%	f.snr_db = inf;
	f.scale_noise = norm(yt) / norm(yn) / exp(f.snr_db / 20);
	yi = yt + f.scale_noise * yn;
	nrms(yi,yt)
	clear xup0 xcp0 yn %xe
%prompt
end


%% "uncorrected" conj. phase recon - no field correction
if ~isvar('xup0'), printm 'xup0'
	xup0 = Gn' * (wi_basis .* yi);
	xup0 = ig.embed(xup0);
	im clf
	im([xtrue; xup0], 'Uncorrected Reconstruction'), cbar
prompt
end


%% slow "exact" conj. phase recon for comparison - with field correctoin
if 0 && ~isvar('xcp0'), printm 'xcp0'
	xcp0 = Ge_zmap' * (wi_basis .* yi);
	xcp0 = ig.embed(xcp0);
	im clf, im([xtrue; xcp0], 'Conjugate Phase Reconstruction'), cbar
prompt
end

% fast conj. phase recon for comparison - with field correctoin
if 0 || ~isvar('xcp1'), printm 'xcp1'
	xcp1 = Gm' * (wi_basis .* yi);
	xcp1 = ig.embed(xcp1);
	im clf, im([xtrue; xcp1], 'Conjugate Phase Reconstruction'), cbar
prompt
end

if 1
	im clf, im([xtrue; xup0; xcp1], 'CP Comparison'), cbar
	nrms(xup0, xtrue) % huge due to spatial distortion along y!
	nrms(xcp1, xtrue)
end


%% penalty
if 0 || ~isvar('R'), printm 'R'
	% scale beta by fov^4 since A'A and 2D.
	f.beta = 2^3;
	R = Reg1(ig.mask, ...
		'type_penal', 'mat', 'type_denom', 'matlab', ... % complex?
		'pot_arg', {'hyper3', 0.1}, 'beta', f.beta);

	if 1 % explore resolution
		[psf var] = qpwls_psf(Gn, R, 1, ig.mask);
		im(psf), cbar
		printm('stddev = %g', sqrt(var * prod(N)))
	prompt
	end
end


%if 0 && ~isvar('Tm'), printm 'Tm: Gmri gram'
%	Tm = build_gram(Gm, 1);
%return
%end

if 0 % check approximation accuracy
% todo: seems not good enough?
	max_percent_diff 'Ge_zmap * xtrue(ig.mask)' 'Gm * xtrue(ig.mask)'
return
end


%% CG
if ~isvar('xcg1'), printm 'xcg1 iterative'
	f.niter = 15;

	xinit = xcp1;
%	xinit = ig.zeros;
	jf_show_iter('mask', ig.mask, 'clim', f.clim); % initialize
	xcg1 = pwls_pcg1(xinit(ig.mask), Gm, 1, yi(:), R, ...
			'userfun', @jf_show_iter, ...
			'niter', f.niter);
	xcg1 = embed(xcg1(:,end), ig.mask);
	im clf, im([xtrue; xcg1], 'xcg1'), cbar

%	if 0 && ~isvar('bb'), printm 'bb'
%		bb = Gm' * yi(:);
%		im clf, im(ig.embed(bb), 'bb'), cbar
%		prompt
%	end
%
%	if 0 && ~isvar('xcg2'), printm 'xcg2'
%		xcg2 = qpwls_pcg2(xinit(ig.mask), Tm, bb, R.C, 'niter', f.niter);
%		xcg2 = ig.embed(xcg2(:,end));
%		im clf, im(xcg2, 'xcg2'), cbar
%		prompt
%	end
prompt
end


if 1
	ix = 10:nx-10;
	im clf
	im([xtrue(ix,:); xup0(ix,:); xcp1(ix,:); xcg1(ix,:)], [0 8]), cbar
	title ''
	nrms(xup0, xtrue) % huge due to spatial distortion along y!
	nrms(xcp1, xtrue)
	nrms(xcg1, xtrue)
	label = {'True', 'Uncorrected', 'Conjugate Phase', 'MBIR'};
	xpos = 2+([1:4]-0.5)*length(ix);
	for ii=1:4
		ir_text(xpos(ii), 0, label{ii}, ...
			'horiz', 'center', 'vertical', 'bottom')
	end
	put = @(k,x) ir_text(xpos(k), ny+1, ...
		sprintf('NRMSE = %.0f\\%%', 100*nrms(x,xtrue)), ...
			'horiz', 'center', 'vertical', 'top');
	put(2, xup0)
	put(3, xcp1)
	put(4, xcg1)
	axis off
%	ir_savefig fig1 % for IEEE SP Magazine paper
return
end

if im % images
	clf, pl=130;
	t = [ [xtrue, xcp1]; [xup0, xcg1] ];
	im('notick', abs(t), f.clim)
	axis off, title ''
	cbar
	tt = @(x,y,s) text(x, y, s, 'horiz', 'center', 'fontsize', 18);
	tt(nx/2, -0.1*ny, 'True');
	tt(3*nx/2, -0.1*ny, sprintf('Uncorrected'));
	tt(nx/2, 2.1*ny, sprintf('Conj. Phase'));
%	tt(nx/2, 2.23*ny, sprintf('L=%d', L));
	tt(3*nx/2, 2.1*ny, sprintf('CG-NUFFT'));
	tt(3*nx/2, 2.23*ny, sprintf('L=%d', L));
end
