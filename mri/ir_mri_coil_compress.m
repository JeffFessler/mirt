 function [odata, sing, Vr] = ir_mri_coil_compress(idata, varargin)
%| MRI coil compression via PCA
%function [odata, S] = ir_mri_coil_compress(idata, varargin)
%|
%| Given multiple MRI surface coil images (idata), use SVD/PCA
%| to find a smaller number of virtual coil images (odata).
%|
%| todo: currently ignores noise correlations
%|
%| in
%|	idata	[(N) n_in]	noisy complex images (2D or 3D) for each coil
%|
%| options
%|	'ncoil'			desired # of virtual coils (default: 1)
%|
%| out
%|	odata	[(N) ncoil]	virtual coil images
%|	sing	[n_in]		singular values
%|	Vr	[n_in, ncoil]	compression matrix for reducing other data
%|
%| Copyright 2016-12-09, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(idata, 'test')
	ir_mri_coil_compress_test
return
end

arg.ncoil = 1;
arg = vararg_pair(arg, varargin);

idim = size(idata);
n_in = idim(end);
idata = reshape(idata, [], n_in); % [*N n_in]
[~, S, V] = svd(idata, 'econ');

Vr = V(:,1:arg.ncoil); % [n_in ncoil] compression matrix with rank = ncoil
odata = idata * Vr; % [*N ncoil] compressed data
odata = reshape(odata, [idim(1:end-1), arg.ncoil]); % [(N) ncoil]

sing = diag(S);


function ir_mri_coil_compress_test

f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
xtrue = single(imread(f.xtrue)');
xtrue = xtrue(2:end-1,2:end-1); % make it 256^2
[nx, ny] = size(xtrue);
%{ 
atrue = 2*pi * (-0.5+([1:nx]'/nx * [1:ny]/ny).^0.5); % smooth phase
xtrue = xtrue .* exp(1i * atrue); % phase
%}

ncoil = 8;
smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 192/nx, ...
	'ncoil', ncoil, 'rcoil', 100);
idata = smap .* repmat(xtrue, [1 1 ncoil]);

snr2sigma = @(db, yb) exp(-db/20) * norm(yb(:)) ...
	/ sqrt(numel(yb)) / sqrt(2); % for complex noise
        sig = snr2sigma(50, idata);
        idata = idata + sig * (randn(size(idata)) + 1i * randn(size(idata)));

nkeep = 4;
[odata, S] = ir_mri_coil_compress(idata, 'ncoil', nkeep);
for iz=1:nkeep
	odata(:,:,iz) = odata(:,:,iz) / max(col(abs(odata(:,:,iz))));
end

if 0 % figures for 551 LR example
	im clf
	im('row', 1, idata, '')
%	cbar
	titlef('MR images from %d coils', 8)
%	ir_savefig cw ir_mri_coil_compress_1a
	im('row', 1, odata, '')
	titlef('First %d singular vectors (virtual coil images)', nkeep)
%	ir_savefig cw ir_mri_coil_compress_1b
	subplot(211)
	style1
	ir_fontsize label 24
	ir_fontsize tick 20
	ir_fontsize title 24
	plot(1:nkeep, S(1:nkeep)/S(1), 'bo', ...
		nkeep+1:ncoil, S(nkeep+1:ncoil)/S(1), 'rx')
	ytick([0 0.2 0.5 1])
	xlabelf('$k$'), ylabelf('$\sigma_k$')
	titlef('Singular values for MRI coil compression example')
	xlim([1-0.2 ncoil+0.2]), ylim([0 1.1]), grid
%	ir_savefig cw ir_mri_coil_compress_1c
return
end

if im
	im plc 3 1
	im('row', 1, smap)
	im('row', 1, idata)
	im pl 3 2
	im(5, 'row', 1, odata)
	im subplot 6
	plot(1:nkeep, S(1:nkeep), 'o', nkeep+1:ncoil, S(nkeep+1:ncoil), 'x')
	titlef('percent kept %.1f\%%', (norm(S(1:nkeep)) / norm(S))^2*100)
end
