 function xpad = ir_pad_into_center(x, npad, varargin)
%function xpad = ir_pad_into_center(x, npad, varargin)
%|
%| Zero pad an input signal x symmetrically around "0" (image center).
%| Useful for DFT/FFT of PSF.
%| using this avoids matlab's padarray.m which is in image proc toolbox
%|
%| option
%|	'circ'	0|1	if 1, circular padding instead of 0.  (default 0)
%| out
%|	xpad	[npad]
%|
%| Originally by A. Yendiki, modified by Jeff Fessler, 2005-7-26.
%| 2015-04-12, renamed, modified to handle case odd npad and for even-sized x

if nargin < 1, ir_usage, end
if nargin == 1 && streq(x, 'test'), ir_pad_into_center_test, return, end

arg.circ = false;
arg = vararg_pair(arg, varargin);

if length(npad) == 1
	if min(size(x)) == 1 % 1d.  kludge: wrong if size(x) = [nx 1 nz]
		if size(x,1) == 1
			npad = [1 npad];
		else
			npad = [npad 1];
		end
	else % n-dimensional; pad all dimensions the same amount
		npad = repmat(npad, [1 ndims(x)]);
	end
end


ndim = ndims(x);
if ndim ~= length(npad)
	error 'Incorrect number of dimensions'
end

if any(size(x) > npad)
	pr size(x)
	pr npad
	fail('npad too small')
end
shift = ceil((npad - size(x))/2);

if arg.circ % circular boundary condition padding
	if any(round(shift) ~= shift)
		fail 'circular needs same padding on both sides'
	end
	args = cell(ndim,1);
	for id = 1:ndim
		nold = size(x,id);
		nnew = npad(id);
		args{id} = 1 + mod([0:nnew-1]-shift(id), size(x,id));
	end
	args = ndgrid_jf('cell', args{:});
	index = sub2ind(size(x), args{:});
	xpad = x(index);

else
	args = cell(ndim,1);
	for id = 1:ndim
		nold = size(x,id);
		nnew = npad(id);
		if nold > nnew
			fail('Padding[%d]=%d too small cf %d', id, nnew, nold)
		end
		offset = ceil((nnew - nold)/2);
		if rem(nnew,2) && ~rem(nold,2) % nnew odd and nold even
			offset = offset - 1; % 2015-04-12 fix
		end
%		pr '[nold nnew offset]'
		args{id} = [1:nold] + offset;
	end
	if islogical(x)
		xpad = false(npad);
	else
		xpad = zeros(npad, class(x));
	end
	xpad(args{:}) = x;
end


function ir_pad_into_center_test

if 1
	jf_equal(ifftshift(ir_pad_into_center([1 2 1], [5])), [2 1 0 0 1])
	jf_equal(ifftshift(ir_pad_into_center([1 2 1], [6])), [2 1 0 0 0 1])
	jf_equal(ifftshift(ir_pad_into_center([1 2 1]', [5]))', [2 1 0 0 1])
	jf_equal(ifftshift(ir_pad_into_center([1 2 1]', [6]))', [2 1 0 0 0 1])
end

% 1d tests for odd and even sizes
for nx=6:7 % odd and even sizes
	for px=2:3 % odd and even sizes
		x2 = 1:nx;
		psf2 = 1:px;
		pad2 = ir_pad_into_center(psf2, size(x2));
		y = convn(x2, psf2, 'same');
		z = round(ifft(fft(x2) .* fft(ifftshift(pad2))));
		jf_equal(y(2:end-1), z(2:end-1))
	end
end

t1 = ir_pad_into_center(psf2, [6]);
t2 = ir_pad_into_center([1 2 1], [7 5]);
t3 = ir_pad_into_center([1 2 1]', [7 3]);
t4 = ir_pad_into_center(ones(3), [5 7]);

x = reshape(1:3*4, 3, 4);
y = ir_pad_into_center(x, [3+4 4+2], 'circ', 1);
if exist('padarray', 'file')
	z = padarray(x, [2 1], 'circular', 'both');
	jf_equal(y, z)
end

if im
	pr t1
	pr t2
	pr t3
	pr t4
	pr x
	pr y
	pr z
end
