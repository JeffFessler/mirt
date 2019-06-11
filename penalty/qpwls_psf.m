 function [psf, var, fwhm, mtf, con] = qpwls_psf(A, C, beta, mask, W, varargin)
%function [psf, var, fwhm, mtf, con] = qpwls_psf(A, C, beta, mask, W, [option])
%|
%| Compute FFT-approximate PSF
%| for quadratically-penalized shift-invariant unweighted least-squares
%| uses psf = [A'WA + beta C'C]^{-1} A'W y, where usually y = A e
%|
%| PSF is computed for pixel at position "offset" relative to
%| the "center" pixel associated with an FFT psf(nx/2+1, ny/2+1)
%| in
%|	C	can be C matrix, or C'Cej, or R object
%|
%| options
%|	'offset'	coordinates of pixel of interest (default in 2d: [0 0])
%|	'dx'		pixel size, to report fwhm in physical units (def: 1)
%|	'dz'		slice spacing, report z.fwhm in physical units (def: 1)
%|	'yb'		substitute this data for A e (default: uses A*e)
%|			this is useful if recon model differs from "truth"!
%|	'shift0'	PSF centered at 0, even if 'offset' (default: 0)
%|	'real_b'	make spectrum of A'*W*yb real, nonnegative,
%|			only applies if yb provided. (default: 0)
%|	'psftype'	'reale' (default) or 'abs' specifies which function
%|			is used to convert psf from complex-valued to real-valued
%|			if psf spectrum is not hermitian symmetric
%|	'fwhmtype'	'slice' (default) or 'profile' (new) or 'fwhm3' (not)
%|			or '2z' to return fwhm of each plane of "3d" psf
%|			or 'none' to prevent fwhm computation
%|	'gtwgej'	user-provided A'*W*A*e_j
%|	'loop'	0|1	if 1, then loop to prompt for other betas.  (default 0)
%|	'cgrad_delta'	default: 1e-2 (todo: document)
%|
%| out
%|	psf
%|	var	approximate relative variance for the center pixel
%|	fwhm	FWHM of psf
%|	mtf	structure with MTF
%|	con	approximate condition number
%|
%| caution: fwhm likely inaccurate if the psf is centered near the image edge
%| yet is wide enough to overlap the edge.  'shift0' may help in such cases.
%|
%| Copyright 2002-1-30, Jeff Fessler, University of Michigan
%| 2015-08-13, changes by Valur Olafsson for FWHM

% todo: option for separating "x" and "y" spatial resolutions in 2d

if nargin == 1 && streq(A, 'test'), qpwls_psf_test; return, end
if nargin < 4, ir_usage, end

arg.chat = 0;
arg.dx = 1; % 1 pixel distance unit
arg.dz = 1; % 1 pixel distance unit
arg.iwarn = 'disp';
arg.offset = zeros(size(size(mask)));
arg.yb = [];
arg.gtwgej = [];
arg.real_b = false;
arg.shift0 = false;
arg.psftype = 'reale';
arg.fwhmtype = 'slice';
arg.loop = false;
arg.cgrad_delta = 1e-2;

% support old-style argument
if length(varargin) && isnumeric(varargin{1})
	arg.offset = varargin{1};
	varargin = {varargin{2:end}};
end

arg = vararg_pair(arg, varargin);

if ~isvar('W') || isempty(W)
	W = 1;	% cheap identity
end

if length(arg.offset) ~= ndims(mask), error 'bad offset dim', end

if ~islogical(mask), warn 'not logical mask', mask = logical(mask); end

%
% unit vector at a pixel
% relative to 'center' pixel in FFT sense
%
dim = size(mask);
ej = zeros(dim);
ej( subv2ind(dim, floor(dim/2) + 1 + arg.offset) ) = 1;
if ~sum(ej(:) .* mask(:))
	warn 'offset outside mask'
	psf = nan(size(mask));
	var = NaN;
	fwhm = NaN;
return
end

% jth column of gram matrix: F ej = A'WA ej
if isempty(arg.gtwgej)
	% make sure power-of-2 dimensions?  nah - make user do it.
	if 0
		dim = size(gtg);
		dpad = 2 .^ ceil(log2(dim-0.1))
		if any(dpad ~= dim)
			warn(['non power of 2, padding ', num2str(dpad)])
		end
	end

	% jth column of Fisher information: A'WA ej
	if issparse(A)
		gtg = embed(A' * (W * (A * sparse(ej(mask(:))))), mask);
	else
		gtg = embed(A' * (W * (A * ej(mask(:)))), mask);
	end

else % user provided (to save time for repeated calls)
	gtg = arg.gtwgej;
end

%
% jth column of regularization matrix: R ej = C'C ej
%
ctc = form_ctc(C, beta, ej, mask, size(gtg), arg.cgrad_delta);

%
% phase due to position offset
%
ndim = length(dim);
if ndim == 2
	k1 = [0:dim(1)-1] / dim(1);
	k2 = [0:dim(2)-1] / dim(2);
	[k1, k2] = ndgrid(k1, k2);
	phi = 2*pi * (k1*arg.offset(1) + k2*arg.offset(2));
elseif ndim == 3
	k1 = [0:dim(1)-1] / dim(1);
	k2 = [0:dim(2)-1] / dim(2);
	k3 = [0:dim(3)-1] / dim(3);
	[k1, k2, k3] = ndgrid(k1, k2, k3);
	phi = 2*pi * (k1*arg.offset(1) + k2*arg.offset(2) + k3*arg.offset(3));
else
	if any(arg.offset ~= 0)
		warn 'phase correction not implemented'
	end
	phi = 0;
end
f.phase = exp(1i*phi);	clear k1 k2 k3 phi

% "spectrum" of Fisher information
f.gtg = fftn(ifftshift(gtg));
f.gtg = f.gtg .* f.phase;
f.gtg = reale(f.gtg, arg.iwarn);
if max(f.gtg(:)) <= 0
	warn 'spectrum of Fisher information must have nonzero maximum!'
	keyboard
end
if min(f.gtg(:)) < 0
	printm('gtg, negative FT: %g%%, n=%d', ...
		min(f.gtg(:)) / max(f.gtg(:)) * 100, sum(f.gtg(:) < 0))

	f.gtg = max(f.gtg,0); % ensure FT(A'A) is >= 0 for inverse filter !!!
end

f.ctc = fftn(ifftshift(ctc));
f.ctc = f.ctc .* f.phase;
f.ctc = reale(f.ctc, arg.iwarn);
if min(f.ctc(:)) < 0
	printm('ctc, negative FT: %g%%, n=%d', ...
		min(f.ctc(:)) / max(f.ctc(:)) * 100, sum(f.ctc(:) < 0))
	f.ctc = max(f.ctc,0); % set negatives to zero
end

if ~isempty(arg.yb)
	gtwy = embed(A' * (W * arg.yb(:)), mask);
	f.gtwy = fftn(ifftshift(gtwy));
	f.gtwy = f.gtwy .* f.phase;
	if arg.real_b
		f.gtwy = reale(f.gtwy, arg.iwarn);
		f.gtwy = max(f.gtwy,0);
	end
	if ~arg.shift0
		f.gtwy = f.gtwy .* conj(f.phase);
	end
end

% now loop to allow user to interactively pick betas
beta_original = beta;
repeater = true;
beta_list = [];
fwhm_list = [];

while repeater
	% hessian "spectrum"
	f.hes = f.gtg + beta * f.ctc;

	con = max(abs(f.hes(:))) / min(abs(f.hes(:)));
	printm('approximate condition number %g', con)

	if any(f.hes(:) <= 0)
		warn('%d zeros in Hessian - results may be meaningless', ...
			sum(f.hes(:) == 0))
		f.hes(f.hes == 0) = inf;
	end

	if ~isempty(arg.yb)
		f.psf = f.gtwy ./ f.hes;
		psf = fftshift(ifftn(f.psf));

	else

		if arg.shift0
			f.psf = f.gtg ./ f.hes;
		else
			f.psf = conj(f.phase) .* f.gtg ./ f.hes;
		end
		psf = fftshift(ifftn(f.psf));
	end

	% note: psf could be complex if spectrum is not hermitian symmetric
	switch arg.psftype
	case 'reale'
		psf = reale(psf, arg.iwarn);
	case 'abs'
		psf = abs(psf);
	otherwise
		fail('unknown psftype %s', arg.psftype)
	end

	if arg.dx == 1
		units = '[pixel]';
	else
		units = '[dis]';
	end

	if streq(arg.fwhmtype, 'none')
		fwhm = [];
	else
		fwhm = qpwls_psf_fwhm(psf, ndim, arg, units); % report fwhm
	end

	% predicted relative variance for white noise based on local Fourier analysis
	var = mean( abs(f.gtg(:)) ./ abs(f.hes(:)).^2 );

	% here is how to compute the local autocovariance function,
	% from the local power spectral density,
	% to within some scale factors (not tested):
	% cov = reale(fftshift(ifftn( f.gtg(:) ./ abs(f.hes).^2 )));

	if 0
		psf = unpadn(psf, dim);
	end

	% if desired, display images, primarily for debug
	if arg.chat && im
		clf, im pl 3 3
		ix = [-dim(1)/2:dim(1)/2-1] - arg.offset(1);
		iy = [-dim(2)/2:dim(2)/2-1] - arg.offset(2);
	%	im(1, ix, iy, abs(gtg), '|A''WA|'), cbar h
		im(1, ix, iy, real(gtg), 'real(A''WA)'), cbar h
	%	im(2, ix, iy, imag(gtg), 'imag(A''WA)'), cbar h
		im(2, ix, iy, ctc, 'C''C'), cbar h
		im(3, ix, iy, psf, 'psf'), cbar h

		ix = [-dim(1)/2:dim(1)/2-1];
		iy = [-dim(2)/2:dim(2)/2-1];
		im(4, ix, iy, fftshift(f.gtg), 'F(A''WA)'), cbar h
		im(5, ix, iy, fftshift(f.ctc), 'F(C''C)'), cbar h
		im(6, ix, iy, fftshift(f.hes), 'F(Hess)'), cbar h

		im(7, ix, iy, fftshift(f.gtg) == 0, 'F(A''WA) == 0'), cbar h
		im(8, ix, iy, fftshift(f.ctc) == 0, 'F(C''C) == 0'), cbar h
		f.hes(f.hes == inf) = 0;
		im(9, ix, iy, fftshift(f.hes) == 0, 'F(Hess) == 0'), cbar h
	end

	if arg.loop
		beta_list(end+1) = beta;
		fwhm_list(end+1,:) = fwhm(:)';
		[beta_list tmp] = sort(beta_list);
		fwhm_list = fwhm_list(tmp,:);
		pr beta_list
		pr fwhm_list
		if 1 && im
			clamp = @(x,a,b) min(max(x, a), b);
			for id=1:ndims(psf)
				ii{id} = ceil(dim(id)/2+1) + [-10:10] + arg.offset(id);
				ii{id} = ii{id}(ii{id} >= 1 & ii{id} <= dim(id));
			end
			clf, im(psf(ii{:}))
		end
		tmp = input('enter new beta scale [q]: ', 's');
		if isempty(tmp) || streq(tmp, 'q')
			repeater = false;
		else
			beta_scale = eval(tmp);
			pr beta_scale
			beta = beta_original * beta_scale;
		end
	else
		repeater = false;
	end
end

%
% mtf
%
if nargout > 3
	mtf.kx = [-dim(1)/2:dim(1)/2-1] / dim(1);
	mtf.ky = [-dim(2)/2:dim(2)/2-1] / dim(2);
	[kx ky] = ndgrid(mtf.kx, mtf.ky);
	mtf.kk = sqrt(kx.^2 + ky.^2);
	mtf.mtf = fftshift(f.psf);
end


% form_ctc_strum()
% in struct/strum cases
function ctc = form_ctc_strum(R, ej, size_gtg, cgrad_delta, mask)

try
	if isvar('R.pot_arg') % check cgrad_delta vs R.arg.delta
		tmp = R.pot_arg;
		tmp = tmp{1}; % todo: why?
		if numel(tmp) == 2 && ischar(tmp{1}) && isnumeric(tmp{2})
			delta = tmp{2};
			if cgrad_delta > tmp{2} / 100
				warn('delta=%g but cgrad_delta=%g.  BAD!', ...
					delta, cgrad_delta)
			end
		end
	end
	ctc = R.cgrad(R, cgrad_delta * ej(mask(:))) / cgrad_delta; % trick
	ctc = embed(ctc, mask);
catch
	error 'bug: unknown Reg type'
end

% trick: allow R to have a different size than A for "extended support"
if any(size(ctc) ~= size_gtg)
	tmp = ctc;
	ext = (size_gtg - size(ctc)) / 2;
	ctc = zeros(size_gtg);
	ctc((1:size(tmp,1)) + ext(1) , (1:size(tmp,2)) + ext(2)) = tmp;
end


% form_ctc()
% determine C' * C * ej
function ctc = form_ctc(C, beta, ej, mask, size_gtg, cgrad_delta)

dim = size(mask);

if isnumeric(C) % sparse (or full) matrix
	if isequal(size(C), dim)
		ctc = C; % here "C" is impulse response of R for center pixel
	else
		if issparse(C)
			ctc = embed(C' * (C * sparse(ej(mask(:)))), mask);
		else
			ctc = embed(C' * (C * ej(mask(:))), mask);
		end
	end

elseif streq(class(C), 'strum') ... % new Reg1 strum
	|| isstruct(C) % Robject

	if beta ~= 1
		printm 'Using beta ~= 1 for R is probably wrong!|||||'
		warn 'ARE YOU SURE YOU KNOW WHAT YOU ARE DOING???????'
	end

	ctc = form_ctc_strum(C, ej, size_gtg, cgrad_delta, mask);

elseif isa(C, 'fatrix2')
	ctc = C' * (C * ej);

elseif isa(C, 'Fatrix')
	ctc = embed(C' * (C * ej(mask(:))), mask);
else
	error 'unknown regularizer'
end


%
% qpwls_psf_fwhm()
%
function fwhm = qpwls_psf_fwhm(psf, ndim, arg, units)

switch ndim
case 2
	try
		fwhm = fwhm2(psf, 'dx', arg.dx, ...
			'type', 'user', ... % 2015-08-13 Valur Olafsson
			'center', floor(size(psf)/2)+1+arg.offset);
		if arg.dx == 1
			printm('fwhm = %g [pixel]', fwhm)
		else
			printm('fwhm = %g [pixel] or %g [dis]', ...
				fwhm/arg.dx, fwhm)
		end
	catch
		warn 'fwhm2 failed, returning NaN for fwhm'
		fwhm = NaN;
	end

case 3
	if streq(arg.fwhmtype, '2z') % each slice of psf separately
		for ii = 1:size(psf,3)
			fwhm(ii) = qpwls_psf_fwhm(psf(:,:,ii), 2, arg, units);
		end
	return
	end

	if ~arg.shift0 && any(arg.offset) && ~streq(arg.fwhmtype, 'fwhm3')
		warn('3d "fwhm" probably wrong without shift0 option')
	end

	try
		dd = [arg.dx arg.dx arg.dz];
		switch arg.fwhmtype

		case 'none'
			fwhm = [];

		case 'profile'
			% new: use fwhm1 for 1d profiles through psf
			tmp = num2cell(floor(size(psf)/2+1));
			[ix iy iz] = deal(tmp{:});
			fwhm = [
				fwhm1(psf(:,iy,iz), 'dx', dd(1))
				fwhm1(psf(ix,:,iz), 'dx', dd(2))
				fwhm1(psf(ix,iy,:), 'dx', dd(3))
			];

		case 'slice' % old: report separate yz, xz, xy fwhm!
			% old: use fwhm2 for 2d slices through psf
			tmp = @(psf, i1, i2) ...
				fwhm2(squeeze(psf), 'dx', dd(i1), 'dy', dd(i2), ...
				'type', 'user', ...
				'center', floor(size(squeeze(psf))/2)+1+arg.offset([i1 i2]));
			warn 'this is the default, but you probably want the ''profile'' option!'
			fwhm = [
				tmp(psf(floor(end/2)+1,:,:), 2, 3)
				tmp(psf(:,floor(end/2)+1,:), 1, 3)
				tmp(psf(:,:,floor(end/2)+1), 1, 2)
			];

		case 'oldslice'
			fwhm = [
				fwhm2(squeeze(psf(floor(end/2)+1,:,:)))
				fwhm2(squeeze(psf(:,floor(end/2)+1,:)))
				fwhm2(squeeze(psf(:,:,floor(end/2)+1)))];

		case 'newslice' % 2015-08-13 Valur Olafsson
			sizePsf = size(psf);
			fwhm = [
				fwhm2(squeeze(psf(floor(end/2)+1,:,:)), ...
					'type', 'user', ...
					'center', floor(sizePsf([2 3])/2)+1+arg.offset([2 3]))
				fwhm2(squeeze(psf(:,floor(end/2)+1,:)), ...
					'type', 'user', ...
					'center', floor(sizePsf([1 3])/2)+1+arg.offset([1 3]))
				fwhm2(squeeze(psf(:,:,floor(end/2)+1)), ...
					'type', 'user', ...
					'center', floor(sizePsf([1 2])/2)+1+arg.offset([1 2]))];

		case 'fwhm3'
			fwhm = fwhm3(psf);
%			im(psf(:,:,floor(end/2)-1+[-3:3])), cbar

		otherwise
			error 'unknown fwhmtype'
		end
		printm('fwhm = %g %g %g %s', fwhm, units)
	catch
		warn 'fwhm(3d) failed, returning NaN for fwhm'
		fwhm = NaN;
	end

otherwise
	fwhm = NaN;
	warn 'fwhm not implemented'
end


%
% self test; also illustrates use
%
function psf = qpwls_psf_test
ig = image_geom('nx', 64, 'ny', 32, 'dx', 1, 'dy', 1);
tmp = ig.mask; tmp(1) = logical(0); ig.mask = tmp; % stress it with mask
sg = sino_geom('par', 'nb', 80, 'na', 40, 'dr', 1, 'strip_width', 'dr');

%sg.plot(ig), prompt

A = Gtomo2_strip(sg, ig);

%R1 = Reg1(ig.mask, 'offsets', [ig.nx-1]); % single offset for Hugo.  bad!
R1 = Reg1(ig.mask);
C = R1.C;
%C = Csparse('maskbar', ig.mask); % test with sparse C

beta = 2^3;
offset = [10 3];
[psf0 var fwhm mtf con0] = ...
	qpwls_psf(A, C, beta, ig.mask, [], 'offset', offset, 'chat', 1);
printm('approximate stddev = %g', sqrt(var))

if 1 % check matched yb version.
	ej = ig.unitv(ig.nx/2+1+offset(1), ig.ny/2+1+offset(2));
	yb1 = A * ej;
	psf1 = qpwls_psf(A, C, beta, ig.mask, [], 'offset', offset, ...
		'fwhmtype', 'none', ...
		'real_b', 0, ... % set to 1 to match perfectly
		'yb', yb1);
	max_percent_diff(psf0, psf1)
prompt
end

if 1
	ell = [offset+0.5 0.5 0.5 0 1];
	yb2 = ellipse_sino(sg, ell, 'oversample', 3);
%	im clf, im(121, yb1), im(122, yb2)
%	max_percent_diff(yb1, yb2)
	psf2 = qpwls_psf(A, C, beta, ig.mask, [], 'offset', offset, ...
		'loop', 0, ...
		'chat', 0, 'shift0', 0, 'yb', yb2);
	max_percent_diff(psf0, psf2)
	im clf, im(cat(3, psf0, psf2, psf0-psf2))
prompt
end

if 0
	im(mtf.kx, mtf.ky, mtf.mtf, 'MTF'), cbar h
	plot(mtf.kk, mtf.mtf, '.'), title 'MTF', axisx(0, 0.7)
end

if 1 % test "R" version
	R3 = Reg1(ig.mask, 'beta', beta);
	psf3 = qpwls_psf(A, R3, 1., ig.mask, [], 'offset', offset, 'chat', 0);
	equivs(psf0, psf3)
end
