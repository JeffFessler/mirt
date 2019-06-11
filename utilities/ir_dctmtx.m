  function out = ir_dctmtx(siz)
%|function out = ir_dctmtx(siz)
%|
%| 1D or 2D DCT orthonormal matrix
%|
%| Reference:
%| Jain, Fundamentals of Digital Image Processing, p. 150.

if nargin ~= 1, ir_usage, end
if streq(siz, 'test'), ir_dctmtx_test, return, end

ndim = numel(siz);
switch(ndim)
case 1
	N = siz;
	k = [0:N-1]';
	n = [0:N-1];
	out = sqrt(2 / N) * cos((pi / 2 / N * k) * (2*n + 1));
	out(1,:) = out(1,:) / sqrt(2);

case 2
	d1 = ir_dctmtx(siz(1));
	d2 = ir_dctmtx(siz(2));
	out = kron(d2, d1);

otherwise
	fail('ndim = %d', ndim)
end


function ir_dctmtx_test
N = 6;
d1 = ir_dctmtx(N);
tmp = d1' * d1;
equivs(tmp, eye(N))
tmp = d1 * d1';
equivs(tmp, eye(N))

M = 4;
d2 = ir_dctmtx([N M]);
tmp = d2' * d2;
equivs(tmp, eye(N*M))
tmp = d2 * d2';
equivs(tmp, eye(N*M))
