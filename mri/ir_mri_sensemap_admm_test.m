%| ir_mri_sensemap_admm_test
%| A test script for the toolbox version of ir_mri_sensemap_admm.
%| Note: This function is not optimized like the one in the paper
%| (uses Fatrix etc.) so the test results will not match those published.
%|
%| See M.J. Allison et al., IEEE TMI, 32(3), 556-564, Mar. 2013.
%|
%| Copyright 2013, Michael Allison, University of Michigan
%| 2015-08 Jeff Fessler, many modifications

%% set parameters
if ~isvar('f')
	f.ns = [0.0275, 0.005]; % noise level (sigma) for [bcoil, lcoil] images
%	f.threshD = 0.13; % allison's original value
	f.threshD = 0.4; % jf value to use only most reliable data
	f.dial = 4; % dilation
	f.ir_imfill1_maskS = true; % apply ir_imfill1 to maskS?
	f.ncoil = 1;
	f.order = 2;
	f.l2b = 5;
	f.niter = 200;
	f.ksADMM = 650;
	f.ku0ADMM = 255;

	if ir_is_octave
		ir_fontsize im_axes 6
		ir_fontsize label 9
		ir_fontsize title 9
	end
end


%% synthetic data (based on M Allison file generateSENSEdata.m)
if ~isvar('ftrue'), printm('ftrue')
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.ftrue = [f.dir 'brainweb_t1.jpg'];
	ftrue = single(imread(f.ftrue)');
	ftrue = ftrue(2:end-1,2:end-1); % make it 256^2
	ftrue = ftrue ./ max(abs(ftrue(:))); % normalize image

	% include synthetic phase component
	[nx ny] = size(ftrue);
	[X, Y] = ndgrid([1:nx] ./ nx, [1:ny] ./ ny);
	ptrue = (Y .* X .* pi) - pi/2;
	ftrue = ftrue .* exp(1i * ptrue);
	clear X Y ptrue
end


%% create synthetic sense maps, noisy body coil and local coil data
if ~isvar('smap'), printm('Simulating sensitivity data')
	smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 192/nx, ...
		'scale', 'ssos_center', ...
		'ncoil', f.ncoil, 'orbit_start', -90, 'rcoil', 100);

	rng(0)
	bcoil = ftrue ...
		+ f.ns(1) * (randn(size(ftrue)) + 1i * randn(size(ftrue)));
	lcoil = smap .* repmat(ftrue, [1 1 f.ncoil]) ...
		+ f.ns(2) * (randn(size(smap)) + 1i * randn(size(smap)));

	% normalize coil images (possibly unnecessary)
	mB = max(abs(bcoil(:)));
	bcoil = bcoil ./ mB;
	lcoil = lcoil ./ mB;
	clear mB
end


% create maskD indicating where body coil data is reliable
% create maskS indicating support for sensitivity map estimation
if ~isvar('maskD'), printm('masks')
	maskD = abs(bcoil) >= f.threshD * max(abs(bcoil(:)));

	if f.dial > 0
		if ir_is_octave
			se = strel('disk', f.dial);
			se = struct(se);
			se = se.nhood;
			maskS = convn(single(maskD), se, 'same') > 0;
		else
			se = strel('disk', f.dial, 0);
			maskS = imdilate(maskD, se);
			maskS = logical(maskS);
		end
		clear se
	else
		maskS = maskD;
	end
	if f.ir_imfill1_maskS
	%	maskS = imfill(maskS, 'holes');
		maskS = ir_imfill1(maskS);
	end
end


if 1 % show data
	slim = minmax(abs(smap) .* maskS)';
	im plc 5 4
	im(1, abs(ftrue), '|ftrue|'), cbar
	im(2, abs(bcoil), '|body|'), cbar
%	im(3, 'colorneg', abs(bcoil) .* maskS, '|body masked|'), cbar
%	im(3, abs(bcoil) .* maskS, '|body masked|'), cbar
	im(3, maskD, 'maskD'), cbar
	im(4, maskS, 'maskS'), cbar
	im(5, abs(smap) .* maskS, '|smap| maskS', slim), cbar
	im(6, abs(lcoil), '|surface|'), cbar
drawnow
end


% RMSE of sensitivity map estimate over support estimation mask
srms_stack = @(ss, smap, mask) ...
	sqrt(mean(abs(ss(mask,:) - repmat(smap(mask), [1 ncol(ss)])).^2));
srms3 = @(shat, smap, mask) srms_stack(reshape(shat, nx*ny, []), smap, mask);
srms2 = @(shat, smap) srms3(shat, smap, maskS);
srms = @(shat) srms2(shat, smap);

%% Compute Cholesky based estimate(s) for this small problem
if ~isvar('schof'), printm('Cholesky estimate with full support')
	[~,~,~,C] = ir_reg_diff_zeroed(size(bcoil), ...
		'order', f.order, 'class', 'sparse');

	bI = col(bcoil .* maskD); % apply data mask (for Cholesky)

	% sparse Cholesky estimate over full square
	printm('-- Computing H');
	H = spdiag(abs(bI.^2), 'nowarn') + 2^(f.l2b) * C' * C;
	rhs = conj(bI) .* lcoil(:);
	printm('-- Computing lhs');
	cpu etic
	schof = H \ double(rhs);
	cpu etoc ': backslash time'
	schof = single(embed(schof, true(size(bcoil))));
end

if 1
	im(9, abs(schof) .* maskS, slim, '$\hat{s}$ cholesky full'), cbar
	im(10, abs(schof - smap) .* maskS, '|err| cholesky full'), cbar
	xlabelf('RMSE %.4f', srms(schof))
drawnow
end


% sparse Cholesky estimate over support (maskS)
if ~isvar('schos'), printm('Cholesky estimate with masked support')
	[~,~,~,C] = ir_reg_diff_zeroed(size(bcoil), 'mask', maskS, ...
		'order', f.order, 'class', 'sparse');
	C = C(:,maskS(:));
	bI = bcoil(maskS) .* maskD(maskS);
	H = spdiag(abs(bI.^2), 'nowarn') + 2^(f.l2b) * C' * C;
	rhs = conj(bI) .* lcoil(maskS);
	schos = H \ double(rhs);
	schos = single(embed(schos, maskS));
	clear H bI rhs C
end

if 1
	im(11, abs(schos), slim, '$\hat{s}$ cholesky mask'), cbar
	im(12, abs(schos - smap) .* maskS, '|err| cholesky mask'), cbar
	xlabelf('RMSE %.4f', srms(schos))
drawnow
end


%% ADMM Estimator for full support
if ~isvar('sadmf'), printm('ADMM method for full support')
	args_admm = {'order', f.order, 'bodycoil', bcoil, 'l2b', ...
		f.l2b, 'niter', f.niter, 'isave', 'all', ...
		'thresh', f.threshD, ...
		'conds', f.ksADMM, 'condu0', f.ku0ADMM, 'maskD', maskD, ...
		'init', 'order3', ...
		'etabtw', 1, 'chat', 1};
	cpu etic
	[sadmf, sinit] = ir_mri_sensemap_admm(lcoil, args_admm{:});
	sinit = reshape(sinit, nx, ny);
	sadmf = reshape(sadmf, nx, ny, []);
	cpu etoc ': admm full'
end

if ~isvar('sadms'), printm('ADMM method for masked support')
	sadms = ir_mri_sensemap_admm(lcoil, args_admm{:}, ...
			'init', sinit .* maskS, 'maskS', maskS);
	sadms = reshape(sadms, nx, ny, []);
end

if 1
	im(7, abs(sinit(:,:,end)) .* maskS, slim, '$\hat{s}$ init'), cbar
	im(8, abs(sinit(:,:,end) - smap) .* maskS, '|err| init'), cbar
	xlabelf('RMSE %.4f', srms(sinit(:,:,end)))

	im(13, abs(sadmf(:,:,end)) .* maskS, slim, '$\hat{s}$ admm full'), cbar
	im(14, abs(sadmf(:,:,end) - smap) .* maskS, '|err| admm full'), cbar
	xlabelf('RMSE %.4f', srms(sadmf(:,:,end)))

	im(15, abs(sadms(:,:,end)) .* maskS, slim, '$\hat{s}$ admm mask'), cbar
	im(16, abs(sadms(:,:,end) - smap) .* maskS, '|err| admm mask'), cbar
	xlabelf('RMSE %.4f', srms(sadms(:,:,end)))
end

%% plot vs truth and non-iterative Cholesky estimator
if im
	im subplot 19
	plot(	0:f.niter, srms(sadmf), 'b-', ...
		0:f.niter, srms(sadms), 'm-'), titlef 'RMSE', legend('full', 'mask')
	axis([0 f.niter 0 0.03]), xtick([0 f.niter]), ytick([0 0.016 0.024])
	im subplot 20
	plot(	0:f.niter, srms2(sadmf, schof), 'b-', ...
		0:f.niter, srms2(sadms, schos), 'm-'), titlef 'RMSD'
	axis([0 f.niter 0 0.03]), xtick([0 f.niter]), ytick([0 0.024])
end
