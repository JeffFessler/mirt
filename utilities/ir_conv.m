 function y = ir_conv(x, psf, varargin)
%function y = ir_conv(x, psf, varargin)
%|
%| ND convolution of signal x with psf.
%| This generalizes convn(x, psf, 'same') by supporting and adjoint
%| and by allowing periodic boundary conditions.
%| For the usual case of an odd-sized psf, the 'center' of the psf is [0,0]
%| so there is no apparently translation of the image.
%|
%| in
%|	x	[(N)]	signal
%|	psf	[(M)]	psf, with ndims(psf) = ndims(x)
%|
%| option
%|	'adj'	0|1	default: 0 (if 1 then adjoint)
%|	'per'	0|1|[]	default: 0 (if 1 then periodic boundary conditions)
%|			if logical array with numel(dims(x)), then selectively
%|			periodic along corresponding directions
%|
%| out
%|	y	[(N)]	x * psf
%|
%| 2012-05-22, Jeff Fessler, Univ. of Michigan
%| 2015-04-12, updated so 'per' version works for all odd/even sizes

if nargin < 1, ir_usage, end
if streq(x, 'test'), ir_conv_test, return, end

arg.per = false;
arg.adj = false;
arg = vararg_pair(arg, varargin);

if arg.adj
	psf = conj(flipdims(psf, 'odd', 1));
end

ndim = ndims(x);

switch arg.per
case false
	y = convn(x, psf, 'same');

case true
	y = ir_conv_per(x, psf, true(1,ndim));

otherwise
	y = ir_conv_per(x, psf, arg.per);
end


% ir_conv_pad()
function x = ir_conv_pad(x, psf, per)
ndim = ndims(x);
if ~islogical(per) || numel(per) ~= ndim
	fail('per must be logical and length %d', ndim)
end

zeross = @(a) zeros(size(a), class(a));

%spad = floor(size(psf) / 2) % pad at start (old)
spad = floor((size(psf) - 1) / 2); % pad at start (2015-04-12)
epad = size(psf) - spad - 1; % pad at end

switch ndim
case 2

	if spad(1) || epad(1)
		pad = spad(1);
		tmp1 = x((end-pad+1):end,:);
		pad = epad(1);
		tmp2 = x(1:pad,:);

		if ~per(1)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end

		x = cat(1, tmp1, x, tmp2);
	end

	if spad(2) || epad(2)
		pad = spad(2);
		tmp1 = x(:,(end-pad+1):end);
		pad = epad(2);
		tmp2 = x(:,1:pad);
		if ~per(2)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(2, tmp1, x, tmp2);
	end

case 3
	if ndims(psf) == 2
		spad = [spad 0];
		epad = [epad 0];
	end

	if spad(1) || epad(1)
		pad = spad(1);
		tmp1 = x((end-pad+1):end,:,:);
		pad = epad(1);
		tmp2 = x(1:pad,:,:);

		if ~per(1)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end

		x = cat(1, tmp1, x, tmp2);
	end

	if spad(2) || epad(2)
		pad = spad(2);
		tmp1 = x(:,(end-pad+1):end,:);
		pad = epad(2);
		tmp2 = x(:,1:pad,:);
		if ~per(2)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(2, tmp1, x, tmp2);
	end

	if spad(3) || epad(3)
		pad = spad(3);
		tmp1 = x(:,:,(end-pad+1):end);
		pad = epad(3);
		tmp2 = x(:,:,1:pad);
		if ~per(3)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(3, tmp1, x, tmp2);
	end

otherwise
	fail('ndim=%d unsupported', ndim)
end


% y = ir_conv_per()
% periodic convolution implemented by circular padding followed by convn
function y = ir_conv_per(x, psf, per)

x = ir_conv_pad(x, psf, per);
y = convn(x, psf, 'valid');


% ir_conv_test()
function ir_conv_test

% all combinations of even/odd for 1D rows and cols
psf1_list = { 1:2, 1:2, 1:3, 1:3, [1:2]', [1:2]', [1:3]', [1:3]' };
x1_list = { 1:6, 1:7, 1:6, 1:7, [1:6]', [1:7]', [1:6]', [1:7]' };

% all combinations of even/odd for 2D
x2_list = cell(2^4,1);
psf2_list = cell(2^4,1);
ii = 1;
for px=2:3
for py=2:3
for nx=7:8
for ny=7:8
	psf2_list{ii} = [1:px]'*[1:py];
	x2_list{ii} = [1:nx]'*[1:ny];
	ii = ii + 1;
end
end
end
end

psf_list = {psf1_list{:}, psf2_list{:}};
x_list = {x1_list{:}, x2_list{:}};

for ii = 1:numel(x_list)
%	pr ii
	psf = psf_list{ii};
	x = x_list{ii};

	% test non-periodic version
	yn = ir_conv(x, psf, 'per', false);
	zn = convn(x, psf, 'same');
	jf_equal(zn, yn)
	if 0
		im plc 3 3
		im(1, yn), cbar, im(2, zn), cbar, im(3, zn-yn), cbar
	end

	% test periodic version
	yp = ir_conv(x, psf, 'per', true);
	psf_centered = ir_pad_into_center(psf, size(x));
	zp = ifftn(fftn(x) .* fftn(ifftshift(psf_centered)));
	if 0
		im(4, yp), cbar, im(5, zp), cbar, im(6, zp-yp), cbar
		im(7, yp - yn), cbar, im(8, zp - zn), cbar
	end
	equivs(zp, yp)
end

if 1 % use Gblur to test the adjoint operation here
	for px=2:3
	for py=2:3
	for nx=6:7
	for ny=6:7
		psf = [1:px]'*[1:py];
		mask = true(nx, ny);
		mask(end) = false; % stress
		A = Gblur(mask, 'psf', psf, 'type', 'conv,per');
		B = Gblur(mask, 'psf', psf, 'type', 'fft,same');
		fA = full(A);
		fB = full(B);
		tA = full(A');
		tB = full(B');
		try
			equivs(fA, fB)
			equivs(tA, tB)
		catch
			ir_conv_show(fA, fB, tA, tB)
			x = zeros(nx, ny);
			x(1) = 1;
			yA = A * x;
			yB = B * x;
			im plc 2 2, im(1, x), im(2, yA), im(3, yB), im(4, yB-yA)
			keyboard
		end
	end
	end
	end
	end

	ir_conv_show(fA, fB, tA, tB)

	test_adjoint(A, 'complex', 1);
end

function ir_conv_show(fA, fB, tA, tB)
im plc 3 3
im(1, fA), cbar
im(2, tA), cbar
im(3, tA - fA', 'tA - fA'''), cbar
im(4, fB), cbar
im(5, tB), cbar
im(6, tB - fB', 'tB - fB'''), cbar
im(7, fA - fB, 'fA - fB'), cbar
im(8, tA - tB, 'tA - tB'), cbar
