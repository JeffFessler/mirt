function out = ir_dct8(x, adj)
%| 2D DCT of each 8x8 block
%| set adj=1 for adjoint (transpose)

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
