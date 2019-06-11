 function [out, r1, ipsf] = Gblob(varargin)
%function [out, r1, ipsf] = Gblob(options)
% generate a KB psf or corresponding Gblur object
%
% options
%	mask		if no mask given, then return psf, not Gblur
%	ndim		2 for 2D
%	J		diameter
%	kb_m		KB order
%	kb_alf		KB alpha (shape)
%	over		over-sampling factor (default 1)
%	fine		fine-sampling factor (default 1)
% out
%	psf	(J+1)^ndim	psf (possibly up-sampled) or Gblur
%	r1	[J*fine+1]	radial sample locations
%	ipsf	(2*J+1)^ndim	inverse psf
%
if nargin == 1 && streq(varargin{1}, 'test'), Gblob_test, return, end
%if nargin < 3, help(mfilename), error(mfilename), end

arg.mask = [];
arg.ndim = 2;
arg.J = 4;
arg.kb_m = 2;
arg.kb_alf = [];
arg.over = 1;
arg.fine = 1;

% options
arg = vararg_pair(arg, varargin);

% table of optimized alpha(ndim,J) values
kb_alf_tab(2,4) = 10.83;

if isempty(arg.kb_alf)
	if arg.ndim <= size(kb_alf_tab,1) ...
		&& round(arg.J) == arg.J ...
		&& arg.J <= size(kb_alf_tab,2) ...
		&& kb_alf_tab(arg.ndim, arg.J) ~= 0 ...
		arg.kb_alf = kb_alf_tab(arg.ndim, arg.J);
	else
		error 'no default kb_alf for that J.  specify!'
	end
end

over = arg.over * arg.fine;
r1 = ([1:arg.over] - (arg.over+1)/2) / arg.over;
r1 = outer_sum([0:arg.J*arg.fine]/arg.fine - arg.J/2, r1)';
r1 = r1(:);
rr = 0;
for id=1:arg.ndim
	rr = outer_sum(rr, r1.^2);
end
rr = sqrt(rr);

psf = kaiser_bessel(rr, arg.J, arg.kb_m, arg.kb_alf);
r1 = downsample2(r1, [arg.over 1]);
if arg.ndim==2
	psf = downsample2(psf, arg.over);
else
	error 'not done'
end
if isempty(arg.mask)
	out = psf;
else
	out = Gblur(arg.mask, 'psf', psf);
end

%
% compute inverse PSF if desired
%
if nargout > 2
	if arg.fine > 1
		psf1 = Gblob(varargin{:}, 'fine', 1);
	else
		psf1 = psf;
	end
	npad = 2^8;
	psf2 = ir_pad_into_center(psf1, npad * ones(1,ndims(psf1)));
	Fpsf = reale(fftshift(fftn(ifftshift(psf2))));
	printf('%s: condition number = %g', mfilename, ...
		max(abs(Fpsf(:))) / min(abs(Fpsf(:))))
	iFpsf = 1 ./ Fpsf; % inverse filter
	ipsf = real(ifftshift(ifftn(fftshift(iFpsf))));
	ndim = ndims(ipsf);
	args = cell(1,ndim);
	ix = npad/2+1+[-arg.J:arg.J];
	for id=1:ndim
		args{id} = ix;
	end
	ipsf = ipsf(args{:});

	tmp = convn(psf, ipsf);
	args = num2cell(floor(size(tmp)/2)+1);
	tmp(args{:}) = tmp(args{:}) - 1;
	printf('%s: inv. psf max error = %g', mfilename, max(abs(tmp(:))))
end

function Gblob_test
%psf = Gblob('over', 1);
%psf = Gblob('over', 2);
fine = 1;
[psf r1 ipsf] = Gblob('over', 1, 'fine', fine);
im(r1, r1, psf)
plot(r1, psf, '.-')
im([ipsf; abs(ipsf) > 1e-5]), cbar
return
dr = r1(2)-r1(1);
fwhm2(psf, dr)
nx = 2^8; ny = nx;
psf2 = ir_pad_into_center(psf, [nx  nx]);
im(psf2)
Fpsf = reale(fftshift(fft2(ifftshift(psf2)))) * dr^2;
im(Fpsf), cbar
contour(Fpsf, [0.5 0.5])
w = [-nx/2:nx/2-1]/nx*2*pi;
plot(w, Fpsf(end/2+1,:), '.'), xaxis_pi '-p -p/2 0 p/4 p/3 p/2 2p/3 p', grid
if fine > 1
	ix = nx/2 + 1 + [-nx/2/fine:nx/2/fine-1];
	Fpsf = Fpsf(ix,ix);
	w = (ix - nx/2-1)/nx*2*pi;
	plot(w, Fpsf(end/2+1,:), '.'), xaxis_pi '-p/3 -p/4 0 p/4 p/3', grid
end
iFpsf = 1 ./ Fpsf; % inverse filter
ipsf = reale(fftshift(ifft2(fftshift(iFpsf))));
ix = -9:9;
ipsf = ipsf(end/2+1+ix,end/2+1+ix);
im(ipsf), cbar
plot(ipsf(ix==0,:))
% im(conv2(ipsf, psf)), cbar
%G = Gblob;
