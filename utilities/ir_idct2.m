 function x = ir_idct2(y)
%function x = ir_idct2(y)
%|
%| Compute 2D inverse DCT of many images
%| (matlab idct2 can handle only one image)
%|
%| in
%|	y [M N (many)]	2D (M by N) DCT coefficients for many images
%|
%| out
%|	x [M N (many)]	corresponding synthesized 2D images
%|
%| 2016 Anish Lahiri
%| 2016-08-024 modified by JF to include test and comments

if nargin < 1, ir_usage, end
if streq(y, 'test'), ir_idct2_test, return, end

[nx, ny, np] = size(y);
x = reshape(y, nx, []); % [M *]
x = idct(x); % [M *] 1D along 1st dim
x = permute(reshape(x, nx, ny, np), [2,1,3]); % [M N np]
x = reshape(x, ny, []); % [N *]
x = idct(x); % [N *] 1D along 2nd dim
x = permute(reshape(x, ny, nx, np), [2 1 3]); % [M N np]
x = reshape(x, size(y));


function ir_idct2_test
m = 8;
n = 10;
nrep = 100;
rng(0);
y = rand(m, n, nrep);
x1 = ir_idct2(y);
x2 = zeros(size(x1));
for ii=1:nrep
	x2(:,:,ii) = idct2(y(:,:,ii));
end
equivs(x1, x2)

y1 = ir_dct2(x1);
equivs(y, y1)

if im % show all the 2D basis functions
	n = 8;
	y = eye(n^2);
	y2 = reshape(y, n, n, []);
	x2 = ir_idct2(y2);
	im plc 1 2
	im(1, x2, '2D: correct')

	x1 = idct(y);
	x12 = reshape(x1, n, n, []);
	im(2, x12, '1D: incorrect')
end
