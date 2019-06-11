  function kern = dsft_gram_hermitify(kern, show)
%|function kern = dsft_gram_hermitify(kern, show)
%|
%| When A is a Gdsft object and W is a Hermitian positive semidefinite matrix,
%| e.g., when W is diagonal with real, nonnegative diagonal elements,
%| the corresponding gram matrix A'WA is positive semidefinite and Toeplitz,
%| and the corresponding kernel h[n] is Hermitian symmetric,
%| i.e., h[-n] = conj(h[n]).
%|
%| This routine takes an input kernel that might be not quite Hermitian
%| and uses averaging to force it to be Hermitian.
%| It is needed because of NUFFT approximation error.
%|
%| in
%|	kern	[(N)]		kernel of gram operator
%|	show	0|1		make plots or images to show kernel?
%|
%| out
%|	kern	[(N)]		kernel of gram operator made to be Hermitian
%|
%| Copyright 2012-06-06, Jeff Fessler, University of Michigan

if nargin == 1 && streq(kern, 'test'), Gdsft_gram test2, clear kern, return, end
if nargin ~= 2, help(mfilename), error(mfilename), end

dim = size(kern);
if numel(dim) == 2 && dim(2) == 1
	dim = dim(1); % trick for 1d
end

switch numel(dim)
case 1
	N1 = size(kern,1);
	n1 = 0:(N1-1);
	m1 = mod(-n1, N1); % circular symmetry
	tmp1 = kern(1+n1);
	tmp2 = kern(1+m1);
	tmp2 = conj(tmp2);
	kern = dsft_gram_hermitify_avg(tmp1, tmp2); % force it to be Hermitian
	dsft_gram_hermitify_show(tmp1, tmp2, show);

case 2
	N1 = size(kern,1);
	N2 = size(kern,2);
	n1 = 0:(N1-1);
	n2 = 0:(N2-1);
	[n1 n2] = ndgrid(n1, n2);
	tmp1 = kern(sub2ind([N1 N2], 1+n1, 1+n2));
	m1 = mod(-n1, N1); % circular symmetry
	m2 = mod(-n2, N2);
	tmp2 = kern(sub2ind([N1 N2], 1+m1, 1+m2));
	tmp2 = conj(tmp2);
	kern = dsft_gram_hermitify_avg(tmp1, tmp2); % force it to be Hermitian
	dsft_gram_hermitify_show(tmp1, tmp2, show);

case 3
	N1 = size(kern,1);
	N2 = size(kern,2);
	N3 = size(kern,3);
	n1 = 0:(N1-1);
	n2 = 0:(N2-1);
	n3 = 0:(N3-1);
	[n1 n2 n3] = ndgrid(n1, n2, n3);
	tmp1 = kern(sub2ind([N1 N2 N3], 1+n1, 1+n2, 1+n3));
	m1 = mod(-n1, N1); % circular symmetry
	m2 = mod(-n2, N2);
	m3 = mod(-n3, N3);
	tmp2 = kern(sub2ind([N1 N2 N3], 1+m1, 1+m2, 1+m3));
	tmp2 = conj(tmp2);
	kern = dsft_gram_hermitify_avg(tmp1, tmp2); % force it to be Hermitian
	dsft_gram_hermitify_show(tmp1, tmp2, show);

otherwise
	fail('dim %d not done', numel(dim))
end


% dsft_gram_hermitify_avg()
function out = dsft_gram_hermitify_avg(in1, in2)
out = (in1 + in2) / 2;
printm('biggest relative change to make Hermitian: %g', ...
	max(col(abs(out - in1))) / max(col(abs(in1))))


% dsft_gram_hermitify_plot1()
function dsft_gram_hermitify_plot1(i, data)
im('subplot', i)
n = numel(data);
if im
	plot([-n/2:n/2-1], fftshift(data), '-o')
	xtick([-n/2 0 n/2-1])
end


% dsft_gram_hermitify_show()
function dsft_gram_hermitify_show(tmp1, tmp2, show)
if ~show, return, end

if size(tmp1,2) == 1 % 1d
	fun = @(n, d) dsft_gram_hermitify_plot1(n, d);
else
	fun = @(n, d) im(n, fftshift(d));
end
err = tmp2 - tmp1; % should be 0
im plc 3 3
fun(1, real(tmp1))
fun(2, real(tmp2))
fun(3, real(err))
fun(4, imag(tmp1))
fun(5, imag(tmp2))
fun(6, imag(err))
fun(7, abs(tmp1))
fun(8, abs(tmp2))
fun(9, abs(err))
prompt
