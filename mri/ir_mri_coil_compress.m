 function [odata, sing, Vr] = ir_mri_coil_compress(idata, varargin)
%| MRI coil compression via PCA
%function [odata, sing, Vr] = ir_mri_coil_compress(idata, varargin)
%|
%| Given data (idata) from multiple MRI surface coils, use SVD/PCA
%| to find data (odata) for a smaller number of virtual coils.
%|
%| todo: currently ignores noise correlations
%|
%| in
%|	idata	[(N) n_in]	noisy complex data for each coil
%|
%| options
%|	'ncoil'		desired # of virtual coils (default: 1)
%|  'thresh'    minimum compression threshold, 0-1 (default: unused)
%|
%| out
%|	odata	[(N) ncoil]	virtual coil data
%|	sing	[n_in]		singular values
%|	Vr	[n_in, ncoil]	compression matrix for reducing other data
%|
%| The method is similar to that of huang:08:asc
%| "A software channel compression technique for faster reconstruction with many channels"
%| MRI Jan. 2008, doi 10.1016/j.mri.2007.04.010
%| But this code does not subtract the data mean.  (Perhaps it should.)
%|
%| Copyright 2016-12-09, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(idata, 'test'), ir_mri_coil_compress_test, return, end

arg.ncoil = 1;
arg.thresh = [];
arg = vararg_pair(arg, varargin);

idim = size(idata);
n_in = idim(end);
idata = reshape(idata, [], n_in); % [*N n_in]
[~, sing, V] = svd(idata, 'econ');

% select number of coils if cutoff set
if arg.thresh
    cum_sing = cumsum(diag(sing)) ./ sum(diag(sing));
    over_ind = find(cum_sing > arg.thresh);
    arg.ncoil = over_ind(1);
end

Vr = V(:,1:arg.ncoil); % [n_in ncoil] compression matrix with rank = ncoil
odata = idata * Vr; % [*N ncoil] compressed data
odata = reshape(odata, [idim(1:end-1), arg.ncoil]); % [(N) ncoil]


% This test case uses *image* data but in practice coil compression
% is applied to k-space data *prior* to image reconstruction!
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

snr2sigma = @(db, yb) 10^(-db/20) * norm(yb(:)) ...
	/ sqrt(numel(yb)) / sqrt(2); % for complex noise
        sig = snr2sigma(50, idata);
        idata = idata + sig * (randn(size(idata)) + 1i * randn(size(idata)));

nkeep = 4;
[odata, S] = ir_mri_coil_compress(idata, 'ncoil', nkeep);
for iz = 1:nkeep
	odata(:,:,iz) = odata(:,:,iz) / max(col(abs(odata(:,:,iz))));
end

if false % figures for 551 LR example
	im clf
	ir_fontsize tick 12
	im('row', 1, idata, '')
%	cbar
	titlef('MR images from %d coils', 8)
	if false % for book
		% ir_savefig cw ir_mri_coil_compress_1a
		xticks([]), yticks([]), title('')
		exportgraphics(gca, 'ir_mri_coil_compress_1a.pdf')
	end
	prompt
	im('row', 1, odata, '')
	titlef('First %d singular vectors (virtual coil images)', nkeep)
	if false % for book
		% ir_savefig cw ir_mri_coil_compress_1b
		xticks([]), yticks([]), title('')
		exportgraphics(gca, 'ir_mri_coil_compress_1b.pdf')
	end
	prompt
	subplot(211)
	style1
	ir_fontsize label 24
	ir_fontsize tick 20
	ir_fontsize title 24
	S = diag(S);
	plot(1:ncoil, S/S(1), 'k')
	hold on
	scatter(1:nkeep, S(1:nkeep)/S(1), 100, 'bo', 'filled')
	scatter(nkeep+1:ncoil, S(nkeep+1:ncoil)/S(1), 100, 'rs', 'filled')
	hold off
	ytick([0 0.2 0.5 1])
	xlabelf('$k$'), ylabelf('$\sigma_k$')
	titlef('Singular values for MRI coil compression example')
	xlim([1-0.2 ncoil+0.2]), ylim([0 1.1]), grid
%	ir_savefig cw ir_mri_coil_compress_1c
%	exportgraphics(gca, 'ir_mri_coil_compress_1c.pdf')
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
	titlef('percent kept %.1f%%', (norm(S(1:nkeep)) / norm(S))^2*100)
end
