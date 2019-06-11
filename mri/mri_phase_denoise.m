 function [x, mask0] = mri_phase_denoise(yi, varargin)
%function [x, mask0] = mri_phase_denoise(yi, [options]) <- recommended usage
%|function [x] = mri_phase_denoise(yi, l2b, niter, chat, wthresh) <- old way
%|
%| in
%|	yi	[(N)]	noisy complex image: y = mag .* exp(1i * x)
%|			can be any size: 1D, 2D, 3D, ...
%|
%| options
%|	l2b		log_2(beta), regularization parameter
%|	order		regularization order (default: 1, for historical
%|				reasons, but 2 is probably preferable)
%|	niter		# of iterations
%|	chat		1 to show pictures
%|	wthresh		fraction of magnitude maximum to include in fitting
%|	init		initial image for iterations
%|	isave		which iterations to save.  (default: last)
%|	'pl'	1|0	1 for new PL method (recommend), 0 for old way (default)
%| out
%|	x	[(N)]	cleaned up phase estimate (radians)
%|	mask0	[(N)]	logical: 1 for high-magnitude pixels
%|
%| Example of weighted "denoising" of MRI phase images.
%| This is a "simple" way to estimate good field inhomogeneity maps
%| from the usual approach of two readouts with a short delay.
%| It also smoothly interpolates over regions with signal voids.
%|
%| Caution: the sign of the field map estimated here is the opposite (negative)
%| of the sign of the field map needed for input to the Gmri object.
%|
%| Copyright 1999, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if ischar(yi)
	[x mask0] = mri_phase_denoise_test(yi, varargin{:});
	if ~nargout, clear x mask0, end
return
end

% defaults
arg.chat = 0;
arg.l2b = -5;
arg.order = 1;
arg.niter = 150;
arg.isave = [];
arg.wthresh = 0.4;
arg.fmax = 0.05; % fraction of max threshold for trimmed median for wi_ml
arg.init = [];
arg.clim = []; % limits for phase display
arg.pl = false; % PL
arg.wi_ml = false; % use wi based on ML instead of threshold

% backword compatible for old argument list: l2b, niter, chat, wthresh
if length(varargin) && isnumeric(varargin{1})
	arg.l2b = varargin{1};
	if length(varargin) >= 2, arg.niter = varargin{2}; end
	if length(varargin) >= 3, arg.chat = varargin{3}; end
	if length(varargin) >= 4, arg.wthresh = varargin{4}; end
else
	arg = vararg_pair(arg, varargin);
end
if isempty(arg.isave), arg.isave = arg.niter; end

mag = abs(yi);
yi = angle(yi);
if arg.chat
	im plc 2 3
	im(1, mag, 'magnitude'), cbar
	im(2, yi, 'raw phase map', arg.clim), cbar
end

dim_yi = size(yi);

%
% specify weights: this needs more work to be automatic!
%
mask0 = mag > arg.wthresh * max(mag(:)); % ignore pixels with "too small" magnitude
if arg.chat
	im(3, mask0, 'weights'), cbar
	im(4, mask0 .* yi, 'masked phase', arg.clim), cbar
end

%mean(mag(:))
%median(mag(:))
%clf, hist(mag(:), 100), pause
%median(mag(mag(:) > 0.05 * max(mag(:))))

if arg.wi_ml
	wi = mri_phase_wi_ml(mag, arg.fmax);
else
	wi = mask0;
end
%W = diag_sp(wi(:));
W = Gdiag(wi);

%
% initial phase image
%
if isempty(arg.init)
	arg.init = yi;
	arg.init(mask0 == 0) = mean(yi(mask0 == 0));
end
if arg.chat
	im(5, arg.init, 'Initial phase', arg.clim), cbar
end

%G = diag_sp(ones(prod(dim_yi),1));
G = Gdiag(ones([dim_yi 1]));

%
% regularizer
%
if arg.order ~= 2, warn('order=2 recommended'), end
mask1 = true(size(yi)); % estimate / extrapolate to *all* pixels
R = Reg1(mask1, 'beta', 2^arg.l2b, 'order', arg.order);
%R = Robject(mask1, 'beta', 2^arg.l2b, 'order', arg.order);
%	'type_denom', 'matlab', ...

if 0 % old way
	[C wjk] = C2sparse('tight', mask1, 8); % todo: cut?
	C = spdiag(sqrt(wjk), 'nowarn') * C; % caution: missing prior to 2005-11-28
	C = sqrt(2^arg.l2b) * C;
end

% report expected blur (at image center)
if 1
	qpwls_psf(G, R, 1, mask1, W);
end

%
% run qpwls algorithm for regularized fitting
%
xinit = arg.init(mask1);
if arg.pl
	med = median(mag(mag(:) > 0.05 * max(mag(:))));
	data = {yi(:), (mag(:)/med).^2}; handle = @phase_dercurv; % PL
%profile on % todo!
	x = pl_pcg_qs_ls(xinit, G, data, handle, R, ...
		'niter', arg.niter, 'isave', arg.isave);
%profile report
	x = embed(x, mask1);
else
	warn 'recommend using "pl" option.  use wls for historical only'
%	data = {yi(:), wi(:)}; handle = @wls_dercurv; % qpwls
	x = qpwls_pcg1(xinit, G, W, yi(:), R.C, ...
		'niter', arg.niter, 'isave', arg.isave);
	x = embed(x, mask1);
end

if 0 % old way
	x = qpwls_pcg(x, G, W, yi(:), 0, R.C, 1, arg.niter);
	x = reshape(x, [dim_yi arg.niter]);
	x = x(:,:,end);
end

if arg.chat
	if ndims(yi) == 2
		im(6, x(:,:,end), 'QPWLS-CG phase', arg.clim), cbar
	else % 3d
		im(6, x(:,:,:,end), 'QPWLS-CG phase', arg.clim), cbar
	end
end


%
% mri_phase_wi_ml()
% wi based on ML estimation
% trick: normalize by median of non-background so that beta is "universal"
%
function wi = mri_phase_wi_ml(mag, fmax)
med = median(mag(mag(:) > fmax * max(mag(:))));
wi = (mag / med).^2;


%
% phase_dercurv()
% wi * (1 - cos(yi - li))
%
function [deriv, curv] = phase_dercurv(data, li, varargin)
yi = data{1};
wi = data{2};
deriv = wi .* sin(li - yi);
curv = wi;


%
% built-in test/example
%
function [xq, mask0] = mri_phase_denoise_test(type, varargin)

% read data
f.dir = path_find_dir('mri');
f.dir = [f.dir '/phase-data/'];
f.dat = [f.dir 'phfit.mat'];
if ~exist(f.dat, 'file')
	fail('edit the path in %s!', mfilename)
end
yi = ir_read_mat(f.dat);
clim = [-0.5 1.5];

%if nargout
%	order = 1;
%else
%	order = 2;
%end
[xq mask0] = mri_phase_denoise(yi, ...
	'clim', clim, 'init', [], 'chat', 1, varargin{:});
%	'init', 5*randn(size(yi)));
%	'init', 5*ones(size(yi)));
xq = xq(:,:,end);
if im
	title 'QPWLS-CG phase (simple wi)'
end

cpu etic
xpl = mri_phase_denoise(yi, 'pl', 1, varargin{:});
cpu etoc 'PL time'
im(4, xpl, 'PL-CG phase', clim), cbar

cpu etic
xml_qpwls = mri_phase_denoise(yi, 'wi_ml', 1, varargin{:});
cpu etoc 'PWLS time'
im(5, xml_qpwls, 'QPWLS-CG (ML wj)', clim), cbar

max_percent_diff(xpl, xml_qpwls)
%max_percent_diff(xpl, xq)
nrms(xml_qpwls, xpl)
nrms(xq, xpl)
%im(4, xpl-xq), cbar

% ir_savefig fig_mr_phase_pl
