function out = ir_dct8(x, adj)
%| 2D DCT of each 8x8 block
%| set adj=1 for adjoint (transpose)

if nargin < 1, ir_usage, end
if nargin < 2, adj = false; end

if streq(x, 'test'), ir_dct8_test, return, end

dc = ir_dctmtx(8);
if adj
	dc = dc';
end

out = zeros(size(x));
[nx ny] = size(x);
for by=1:floor(ny/8)
	for bx=1:floor(nx/8)
		ix = (bx-1)*8+(1:8);
		iy = (by-1)*8+(1:8);
		block = x(ix,iy);
		tmp = (dc * block) * dc';
		out(ix,iy) = tmp;
	end
end

end

function ir_dct8_test
	nx = 32;
	ny = 16;
	X = rand(nx,ny);
	d1 = ir_dct8(X);
%	im(d1)

	b8 = dctmtx(8);
	bx = kron(eye(nx/8), b8);
	by = kron(eye(ny/8), b8);
	d2 = bx * X * by';
%	im(d2)
	assert(max_percent_diff(d1, d2) < 1e-12)

	d3 = kron(by, bx) * X(:);
	d3 = reshape(d3, nx, ny);
%	im(d3)
	assert(max_percent_diff(d1, d3) < 1e-12)
end
