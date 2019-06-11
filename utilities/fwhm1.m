  function [fw, hr, hl] = fwhm1(psfs, varargin)
%|function [fw, hr, hl] = fwhm1(psfs, options)
%|
%| compute fwhm of point-spread function centered at pixel imid
%| and half-right width and half-left width (fw = hr + hl)
%|
%| in
%|	psfs	[np nc]	column psf(s)
%|
%| option
%|	imid	[1]	which pixel is the middle (default: peak)
%|	dx	[1]	pixel size (default: 1)
%|	min0	0|1	if 1, then negative psf value set to 0 (default: 1)
%|	chat	0|1	if 1, then plot
%| out
%|	fw	[nc]	fwhm of each psf column
%|
%| 2012-09-13, added 'min0' option to avoid fwhm < 1

if ~nargin, ir_usage, end
if streq(psfs, 'test'), fwhm1_test, return, end

arg.chat = 0;
arg.dx = 1;
arg.imid = [];
arg.min0 = true;
arg = vararg_pair(arg, varargin);

psfs = squeeze(psfs); % in case [1 1 np]
[np nc] = size(psfs);
if (np == 1) % single row
	psfs = psfs';
	[np nc] = size(psfs);
end

warned = false;
for ic = 1:nc
	psf = psfs(:,ic);
	if isempty(arg.imid)
		imid = imax(psf);
	else
		imid = arg.imid;
	end

	% normalize
	psf = psf / psf(imid);
	if ~warned && (1 ~= max(psf))
		warn 'peak not at center'
		warned = true;
	end

	% right
	ir = sum(cumprod(double(psf((imid+1):np) >= 0.5)));
	if (imid + ir == np)
		hr(ic,1) = ir;
	else
		high = psf(imid + ir);
		ilow = imid + ir + 1;
		low = psf(ilow);
		if arg.min0 && low < 0
			warn('psf(%d) = %g < 0, so using 0', ilow, low)
			low = 0;
		end
		hr(ic,1) = ir + (high - 1/2) / (high-low);
	end

	% left
	il = sum(cumprod(double(psf((imid-1):-1:1) >= 0.5)));
	if (il == imid-1)
		hl(ic,1) = il;
	else
		high = psf(imid - il);
		ilow = imid - il - 1;
		low = psf(ilow);
		if arg.min0 && low < 0
			warn('psf(%d) = %g < 0, so using 0', ilow, low)
			low = 0;
		end
		hl(ic,1) = il + (high - 1/2) / (high-low);
	end
end

hr = hr * arg.dx;
hl = hl * arg.dx;
fw = hr + hl;

if arg.chat && im
	plot(([1:np]-imid)*arg.dx, psf, '-o', ...
	[-hl hr], [0.5 0.5], '-')
	xlabel 'x', ylabel 'psf(x)'
	title(sprintf('fwhm=%g', fw))
end


% fwhm1_test
function fwhm1_test
dx = 3;
nx = 100;
xx = [-nx/2:nx/2-1]' * dx;
fx = 30;
sx = fx / sqrt(log(256));
psf = exp(-((xx/sx).^2)/2);
fw = fwhm1(psf, 'dx', dx, 'chat', 1);

psf = zeros(7,1); psf(5) = 2; psf([4 6]) = -0.2;
fw = fwhm1(psf, 'dx', dx, 'chat', 1);
