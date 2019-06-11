 function [fw, hr, hl] = fwhm(psfs, ii, doplot)
%function [fw, hr, hl] = fwhm(psfs, ii, doplot)
% compute fwhm of point-spread function centered at pixel ii
% and half-right width and half-left width (fw = hr + hl)

if ~nargin, ir_usage, end

warning 'fwhm is obsolete.  use fwhm1 for 1D fwhm'

[np nc] = size(psfs);
if (np == 1)
	psfs = psfs';
	[np, nc] = size(psfs);
end

warned = false;
for ic = 1:nc
	psf = psfs(:,ic);
	if (nargin < 2)
		ii = imax(psf);
	end

	% normalize
	psf = psf / psf(ii);
	if ~warned && (1 ~= max(psf))
		warning('peak not at center')
		warned = true;
	end

	% right
	ir = sum(cumprod(double6(psf((ii+1):np) >= 0.5)));
	high	= psf(ii + ir);
	low	= psf(ii + ir + 1);
	hr(ic,1) = ir + (high - 1/2) / (high-low);

	% left
	il = sum(cumprod(double6(psf((ii-1):-1:1) >= 0.5)));
	high	= psf(ii - il);
	low	= psf(ii - il - 1);
	hl(ic,1) = il + (high - 1/2) / (high-low);
end

fw = hr + hl;

if nargin > 2
	plot(1:np, psf, 'o', ii+[-left right], [0.5 0.5], '-')
end
