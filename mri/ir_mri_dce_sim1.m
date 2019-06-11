%| ir_mri_dce_sim1
%|
%| Generate DCE k-space data for parallel MRI
%|
%| 2014-08-25 Jeff Fessler and Mai Le, University of Michigan

if ~isvar('Ncoil') % parameters
	Ncoil = 8;
	Nframe = 12;	% increase this for better time resolution;
			% keep small for quick debugs
end

if ~isvar('dyn_obj') % dynamic object with 'fine' time sampling as "truth"
	[dyn_obj dce] = ir_mri_dce_obj1('chat', 0);

	im plc 3 3
	im(1, dce.labels, 'labels')

	im subplot 2
	plot(dce.ti, dce.sig, '.-')
	xlabel 'time [min]'
	ylabel 'MR signal'
end


if ~isvar('mask') % spatial mask
	mask = conv2(dce.labels, ones(3), 'same') > 0;
	im(3, mask)
end


if ~isvar('samp1') % k-space dynamic sampling
	[samp1 samp2] = ir_mri_dce_samp1(size(dyn_obj), ...
		'n_tr_merge', dce.n_tr_merge, 'Nframe', Nframe, 'chat', 1);
	tmp = sprintf('samp(1:4:%d)', Nframe);
	im(4, samp1(:,:,1:4:end), tmp)
end


if ~isvar('smap') % coil patterns
	nx = size(dyn_obj,1);
	ny = size(dyn_obj,2);
	smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, ...
		'ncoil', Ncoil, 'coil_distance', 1.5);
	tmp = repmat(mask, [1 1 Ncoil]) .* smap;
	im(5, abs(tmp), '|smap| masked')
end


if ~isvar('kspace_true') % true kspace data [nx ny Nframe Ncoil] (with zeros)
	kspace_true = ir_mri_dce_kspace1(dyn_obj, samp2, ...
			'Nframe', Nframe, 'smap', smap, 'chat', 0);
	pr size(kspace_true)

	ir_show = 1:4:Nframe;
	tmp = log(abs(kspace_true(:,:,ir_show,1))); % 1st coil
	tmp(~samp1(:,:,ir_show)) = -8; % avoid warning about infty
	im(6, tmp, [-8 0] + max(tmp(:)), 'log(|kspace(1:4:?)|)')
end


if ~isvar('xtrue') % 'xtrue' approximation based on frame-averaged
	xtrue = reshape(dyn_obj, [nx ny size(dyn_obj,3)/Nframe Nframe]);
	xtrue = squeeze(mean(xtrue, 3)); % [nx ny Nframe]
	im(7, xtrue)
end


if ~isvar('yi') % noisy kspace data % [Npf Ncoil Nframe] (Npf = # PE per frame)
	ytrue = masker(kspace_true, samp1); % [Ntr Ncoil] (compact - no zeros)
	ytrue = reshape(ytrue, [], Nframe, Ncoil);
	ytrue = permute(ytrue, [1 3 2]); % [Npf Ncoil Nframe]
	if any(ytrue(:) == 0), error 'bug?', end
	snr2sigma = @(db, yb) exp(-db/20) * norm(ytrue(:)) ...
		/ sqrt(numel(ytrue)) / sqrt(2); % for complex noise

	rng(0)
	snr_db = 50;
	sig = snr2sigma(snr_db, ytrue);
	yi = ytrue + sig * (randn(size(ytrue)) + 1i * randn(size(ytrue)));

	tmp = squeeze(yi(:,1,:)); % [Npf Nframe] 1st coil
	tmp = embed(tmp(:), samp1); % [(N) Nframe]
	tmp = log(abs(tmp(:,:,ir_show))); % [(N) frames]
	tmp(~samp1(:,:,ir_show)) = -8; % avoid warning about infty
	im(8, tmp, [-8 0] + max(tmp(:)), 'log(|yi(1:4:?)|)')
drawnow
end


if ~isvar('Acoil') % fatrix for coil maps
	tmp = cell(Ncoil, 1);
	for ic=1:Ncoil
		diag = masker(smap(:,:,ic), mask);
		tmp{ic} = Gdiag(diag, 'mask', mask); % [(N)] x [(N)]
	end
	Acoil = vertcat(tmp{:}); % [(N) Ncoil] x [(N)]
	Acoil = kronI(Nframe, Acoil); % [(N) Ncoil Nframe] x [(N) Nframe]
end

if ~isvar('A') % overall system matrix as fatrix

	% fatrix for time series where each frame may have different sampling
	tmp = cell(Nframe, 1);
	for ir=1:Nframe
		tmp{ir} = Gdft('mask', mask, 'samp', samp1(:,:,ir), ...
			'fftshift', true, 'ifftshift', true); % [Npf] x [(N)]
		tmp{ir} = kronI(Ncoil, tmp{ir}); % [Npf Ncoil] x [(N) Ncoil]
	end
	Aft = block_diag(tmp{:}); % [Npf Ncoil Nframe] x [(N) Ncoil Nframe]

	A = Aft * Acoil; % [Npf Ncoil Nframe] x [(N) Nframe]
end

if 0 % test forward
	tmp = A * xtrue;
	tmp = permute(tmp, [1 2 4 3]);
	equivs(tmp, kspace_true, 'fail', false) % only 1% difference!
return
end

if 0 % test back
	tmp = A' * yi(:);
	tmp = reshape(tmp, sum(mask(:)), Nframe);
	tmp = embed(tmp, mask);
	im(9, abs(tmp))
return
end

if ~isvar('xhat') % run simple LS algorithm to check sizes etc.
	xinit = zeros(size(xtrue), 'single'); % yuch
	tmp = masker(xinit, mask); % convert to masked columns for each frame
	niter = 20;
	xhat = qpwls_pcg1(tmp(:), A, 1, yi(:), 0, 'niter', niter); % column work
	xhat = reshape(xhat, [], Nframe); % convert back
	xhat = embed(xhat, mask); % [(N) Nframe]
	im(9, abs(xhat))
prompt
end

% examine some time curves
% my unregularized LS recon is not very good for rapid enhancing lesion
% perhaps because too few iterations or too poor time resolution
if 1
	labs = unique(dce.labels(:));
%	labs(labs == 0) = []; % not background
	labs = 11:13; % lesions
	rois = zeros(numel(labs), Nframe, 'single');
	for il = 1:numel(labs)
		tmp = dce.labels == labs(il);
		roi = repmat(tmp, [1 1 Nframe]) .* xhat;
		rois(il,:) = sum(sum(roi,1),2) ./ sum(tmp(:)); % mean
	end
	rois = abs(rois); % ?

	clf
	plot(dce.ti, dce.sig(labs,:), 'g.-')
	tt = ([1:Nframe]-0.5) / Nframe * dce.duration_s / 60; % [min]
	hold on
	plot(tt, rois, 'r-o')
	hold off
return
end
