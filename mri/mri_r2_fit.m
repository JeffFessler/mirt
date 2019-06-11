  function r2 = mri_r2_fit(images, telist, varargin)
%|function r2 = mri_r2_fit(images, telist, varargin)
%|
%| fit R2=1/T2 maps to images with different echo times
%| in
%|	images	[(nd) nt]	images with nt different echo times
%|	telist	[nt]		echo times
%|
%| options
%|	'how'		char	'log-ls' - ordinary LS fit to log data
%|	'threshold'	[0,100]	% of maximum signal to keep.  default 0
%|
%| out
%|	r2	[(nd)]		R2 maps
%|
%| 2012-05-16, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(images, 'test'), mri_r2_fit_test, return, end

arg.how = 'log-ls'; % ordinary LS fit to log data
arg.threshold = 0; % what % of maximum signal to keep
arg = vararg_pair(arg, varargin);

dims = size(images);
if numel(telist) ~= dims(end)
	fail('numel(telist) = %d but last dimension of images is %d', ...
		numel(telist), dims(end))
end

switch arg.how
case 'log-ls'
	zlog = log(abs(images));
	zlog = reshapee(zlog, [], numel(telist)); % [*nd nt]

	telist = telist(:); % [nt 1]
	temean = mean(telist);
	tmp = (temean - telist) / sum(telist.^2 - temean^2); % [nt 1]
	r2 = zlog * tmp; % [*nd 1]
	r2 = reshape(r2, [dims(1:end-1) 1]); % [(nd)]
otherwise
	fail('unknown how "%s"', arg.how)
end

if arg.threshold > 0
	tmp = sum(abs(images), numel(dims));
	tmp = tmp / max(tmp(:));
	mask = tmp > (arg.threshold / 100);
	r2 = r2 .* mask;
end


% mri_r2_fit_test
function mri_r2_fit_test

ig = image_geom('nx', 64, 'ny', 60, 'dx', 1);
ell = [0 0 28 28 0 0.05; 11 0 10 10 0 0.05; -10 0 8 8 0 -0.03];
r2true = ellipse_im(ig, ell, 'oversample', 2);
ell = [0 0 28 28 0 1; 0 0 25 25 0 -1; 0 0 23 23 0 1];
pdtrue = ellipse_im(ig, ell, 'oversample', 2);
im plc 3 3
im(1, pdtrue), cbar
im(2, r2true), cbar

telist = [0 10 20 30]; % ms
ztrue = zeros([size(pdtrue) numel(telist)]);
for it=1:numel(telist)
	ztrue(:,:,it) = pdtrue .* exp(-r2true * telist(it));
end

im pl 3 1
im('row', 1, 2, ztrue), cbar

rng(0)
snr = 50;
snr2sigma = @(db, yb) exp(-db/20) * norm(yb(:)) / sqrt(numel(yb)) / sqrt(2); % for complex noise
sigma = snr2sigma(snr, ztrue);
zi = ztrue + sigma * (randn(size(ztrue)) + 1i * randn(size(ztrue)));
printm('snr = %g dB', 20*log(norm(ztrue(:)) / norm(col(ztrue-zi))))

im('row', 1, 3, zi), cbar

r2ls = mri_r2_fit(zi, telist);

clim = minmax(r2true)';
im pl 3 3
im(3, r2ls, clim), cbar
