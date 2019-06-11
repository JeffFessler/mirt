 function ob = Gblur2(mask, varargin)
%function ob = Gblur(mask, options)
%|
%| Construct Gblur object for image restoration.
%|
%| See Gblur_test.m for example usage.
%|
%| in
%|	mask	size(image)	logical array of object support.
%|
%| options
%|	'chat'		verbose printing of debug messages
%|	'psf'		point spread function (aka impulse response)
%|	'type'		type of blur:
%|				'conv,same'	usual case (default)
%|				'conv,per'	(periodic end conditions)
%|				'fft,same'	(periodic end conditions)
%|				todo: allow replicated end conditions!
%|				'imfilter,same'	(requires image toolbox)
%|				'imfilter,circ'	(requires image toolbox)
%|				'imfilter,mirror' (requires image toolbox)
%|					(do not use - adjoint fails)
%|				'sparse'	todo: return sparse matrix
%|	'imfilter_options'	options to 'imfilter' (if used)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = prod(size(mask)) for 'conv,same' type
%|
%| Copyright 2005-4-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mask, 'test'), Gblur_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;

% option defaults
arg.chat = 0;
arg.psf = 1; % identity matrix
arg.type = 'conv,same';
arg.class = 'fatrix2';
arg.imfilter_options = {};

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

if ndims(arg.psf) ~= ndims(mask), error 'psf dim mismatch', end
if any(size(arg.psf) > size(mask)), error 'psf too large', end

if streq(arg.type, 'imfilter,circ') && exist('imfilter') ~= 2
	warn 'no imfilter (no image processing toolbox) so using fft'
	arg.type = 'fft,same';
end

if streq(arg.type, 'imfilter,mirror') && exist('imfilter') ~= 2
	fail 'no imfilter (no image processing toolbox)'
end

if streq(arg.type, 'imfilter,same') && exist('imfilter') ~= 2
	warn 'no imfilter (no image processing toolbox) so using slower convn'
	arg.type = 'conv,same';
end

arg.psf = Gblur_mask_psf_odd(arg.psf);

psf_flip = conj(flipdims2(arg.psf, 'odd', 1));
does_many = false;

switch arg.type

case 'conv,same'
	arg.odim = size(mask);
	arg.psf_flip = psf_flip;
	handle_forw = @(arg,x) convn(x, arg.psf, 'same');
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			convn(y, arg.psf_flip, 'same'));
	does_many = true;

case 'conv,per'
	arg.odim = size(mask);
	arg.psf_flip = psf_flip;
	handle_forw = @(arg,x) ir_conv(x, arg.psf, 'per', true);
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			ir_conv(y, arg.psf_flip, 'per', true));

case 'fft,same'
	arg.odim = size(mask);
	[arg, handle_forw, handle_back] = Gblur_setup_fft_same(arg);

case 'imfilter,circ'
	arg.odim = size(mask);
	arg.psf_flip = psf_flip;
	handle_forw = @(arg,x) imfilter(x, arg.psf, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			imfilter(y, arg.psf_flip, 'conv', 'same', ...
				'circular', arg.imfilter_options{:}));

case 'imfilter,mirror'
	arg.odim = size(mask);
	arg.psf_flip = psf_flip;
	handle_forw = @(arg,x) imfilter(x, arg.psf, 'conv', 'same', ...
		'symmetric', arg.imfilter_options{:});
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			imfilter(y, arg.psf_flip, 'conv', 'same', ...
				'symmetric', arg.imfilter_options{:}));

case 'imfilter,same'
	arg.odim = size(mask);
	arg.psf_flip = psf_flip;
	handle_forw = @(arg,x) imfilter(x, arg.psf, 'conv', 'same', ...
		arg.imfilter_options{:});
	handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			imfilter(y, arg.psf_flip, 'conv', 'same', ...
				arg.imfilter_options{:}));

case 'sparse'
	ob = Gblur_sparse(arg.psf, arg.idim);
	ob = ob(:,arg.mask(:));
	return

otherwise
	error 'unknown blur type'
end

switch arg.class

case 'Fatrix'

	arg.nd = prod(arg.odim);
	arg.np = sum(mask(:));
	dim = [arg.nd arg.np]; % trick: make it masked by default!

	ob = Fatrix(dim, arg, 'caller', 'Gblur', ...
		'abs', @Gblur_abs, ...
		'forw', @Gblur_forw_Fatrix, 'back', @Gblur_back_Fatrix);

case 'fatrix2'

	ob = fatrix2('idim', arg.odim, 'arg', arg, 'does_many', does_many, ...
		'abs', @Gblur_abs, 'forw', handle_forw, 'back', handle_back);

otherwise
	fail('bug')
end


% Gblur_mask_psf_odd(psf)
% having an odd-sized psf facilites adjoint
function psf = Gblur_mask_psf_odd(psf)
nd = ndims(psf);
for id=1:nd
	sz = size(psf, id);
	if ~mod(sz, 2) % even
		warn('size(psf, %d) = %d is even; appending 0', id, sz)
		tmp = size(psf);
		tmp(id) = 1;
		z = zeros(tmp);
		%psf = cat(id, psf, z);
		psf = cat(id, z, psf);
	end
end


% Gblur_sparse()
% a_ij = b[ n(i) - n(j), m(i) - m(j) ]
% n(i) = (i-1) mod N
% m(i) = floor((i-1) / N)
function sp = Gblur_sparse(psf, idim)
if numel(idim) ~= 2
	fail 'not done'
end
fail 'todo: not done'
%sp = sparse(i, j, s, length(ii), size(ob.arg.G, 2));


% Gblur_abs(): |A| for abs(A)
function ob = Gblur_abs(ob)
ob.arg.psf = abs(ob.arg.psf);
if isfield(ob.arg, 'psf_flip')
	ob.arg.psf_flip = abs(ob.arg.psf_flip);
end
if isfield(ob.arg, 'psf_fft') % fixed 2012-04-08
	ob.arg.psf_fft = Gblur_setup_fft_same_fft(ob.arg.psf, ob.mask);
end


% Gblur_setup_fft_same()
function [arg, handle_forw, handle_back] = Gblur_setup_fft_same(arg)

arg.psf_fft = Gblur_setup_fft_same_fft(arg.psf, arg.mask);

handle_forw = @(arg,x) ifftn(fftn(x) .* arg.psf_fft);
handle_back = @(arg,y) fatrix2_maskit(arg.mask, ...
			ifftn(fftn(y) .* conj(arg.psf_fft)));


% Gblur_setup_fft_same_fft()
function psf_fft = Gblur_setup_fft_same_fft(psf, mask)

% % put psf in center of array
% tmp = zeros(size(mask));
% switch ndims(psf)
% case 2
% 	[n1 n2] = size(mask);
% 	[p1 p2] = size(psf);
% %	h1 = (size(psf,1)-1)/2;
% %	h2 = (size(psf,2)-1)/2;
% %	i1 = floor(n1/2) + 1 + [-h1:h1]
% %	i2 = floor(n2/2) + 1 + [-h2:h2]
% 	i1 = floor(n1/2) + 1 + (1:p1) - floor((p1+1)/2);
% 	i2 = floor(n2/2) + 1 + (1:p2) - floor((p2+1)/2);
% 	tmp(i1,i2) = psf;
% case 3
% 	if any(~rem(size(psf),2)), error 'psf size must be odd', end
% 	h1 = (size(psf,1)-1)/2;
% 	h2 = (size(psf,2)-1)/2;
% 	h3 = (size(psf,3)-1)/2;
% 	i1 = floor(size(mask,1)/2) + 1 + [-h1:h1];
% 	i2 = floor(size(mask,2)/2) + 1 + [-h2:h2];
% 	i3 = floor(size(mask,3)/2) + 1 + [-h3:h3];
% 	tmp(i1,i2,i3) = psf;
% otherwise
% 	error 'only 2d & 3d done'
% end
% psf_fft = fftn(fftshift(tmp));

% pad psf with 0, and make size(psf)==size(mask)
sz_mask = size(mask);
sz_psf = ones(1,ndims(mask));
for ii=1:length(sz_psf), sz_psf(ii) = size(psf,ii); end
if any(sz_mask>sz_psf)
    psf = padarray(psf, sz_mask-sz_psf, 'post');
end
center = ceil(sz_psf/2-1/2) + 1; % or center = floor(sz_psf/2) + 1;
psf_fft = fftn( circshift(psf, 1-center) ); % by yusheng


% Gblur_forw_Fatrix(): y = A * x
function y = Gblur_forw_Fatrix(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

switch arg.type
case 'conv,same'
	y = convn(x, arg.psf, 'same');
case 'conv,per'
	y = ir_conv(x, arg.psf, 'per', true);
case 'fft,same'
	if ndims(x) > ndims(arg.mask);
		y = fft_conv_multi(x, arg.psf_fft);
	else
		y = ifftn(fftn(x) .* arg.psf_fft);
	end
case 'imfilter,circ'
	y = imfilter(x, arg.psf, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
case 'imfilter,same'
	y = imfilter(x, arg.psf, 'conv', 'same', arg.imfilter_options{:});
otherwise
	error 'bug'
end

y = ei.shape(y);


% Gblur_back_Fatrix(): x = A' * y
function x = Gblur_back_Fatrix(arg, y)

[y eo] = embed_out(y, arg.odim);

switch arg.type
case 'conv,same'
	x = convn(y, arg.psf_flip, 'same');
case 'conv,per'
	x = ir_conv(y, arg.psf_flip, 'per', true);
case 'fft,same'
	if ndims(y) > ndims(arg.mask);
		x = fft_conv_multi(y, conj(arg.psf_fft));
	else
		x = ifftn(fftn(y) .* conj(arg.psf_fft));
	end
case 'imfilter,circ'
	x = imfilter(y, arg.psf_flip, 'conv', 'same', ...
		'circular', arg.imfilter_options{:});
case 'imfilter,same'
	x = imfilter(y, arg.psf_flip, 'conv', 'same', arg.imfilter_options{:});
otherwise
	error 'bug'
end

x = eo.shape(x, arg.mask, arg.np);


% fft_conv_multi()
function y = fft_conv_multi(x, H)
dim = size(x);
y = zeros(size(x));
y = [];
for ii=1:dim(end)
	tmp = ifftn(fftn(stackpick(x, ii)) .* H);
	y = cat(length(dim), y, tmp); % slow: keeps enlarging y.
end
