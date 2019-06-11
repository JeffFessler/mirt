 function [T, reuse] = Gdsft_gram(ob, W, reuse)
%function [T, reuse] = Gdsft_gram(A, W, reuse)
%|
%| Construct Toeplitz gram matrix object T = A'WA for a Gdsft object.
%| This object performs T*x rapidly using an over-sampled FFT.
%|
%| in
%|	A	[M np]		Gdsft object (fatrix2)
%|	W	[M M]		W = Gdiag() for fatrix2 (often simply "1" or [])
%|				assume W=diag(wi) with real wi so T is Hermitian
%|	reuse	struct		stuff from the last call that can be reused
%|
%| out
%|	T	[np np]		fatrix2 object
%|	reuse	struct		stuff that can be reused on future calls
%|
%| Copyright 2012-06-05, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ob, 'test'), Gdsft_gram_test, return, end
if nargin == 1 && streq(ob, 'test1'), Gdsft_gram_test1, return, end
if nargin == 1 && streq(ob, 'test2'), Gdsft_gram_test2, return, end
if nargin == 1 && streq(ob, 'test3'), Gdsft_gram_test3, return, end
if nargin < 2 || nargin > 3, ir_usage, end

if ~isvar('reuse'), reuse = []; end

if ~isa(ob, 'fatrix2')
	fail('A wrong class: %s', class(ob))
end

[nd np] = size(ob);
forw = @(arg, x) Gdsft_gram_mult(x, arg.fftkern, arg.mask);

% extract wi
if isnumeric(W) % includes case where W is empty
	wi = W;
elseif isa(W, 'fatrix2') && streq(W.caller, 'Gdiag')
	wi = W.arg.diag;
else
	error 'W must be diag_sp or Gdiag or wi array'
end
clear W

if isempty(wi)
	wi = ones(nd, 1); % the usual unweighted case
end

if isscalar(wi)
	wi = wi * ones(nd, 1);
end

if isreal(wi)
	back = forw; % trick: because Hermitian!
else
	fail('not implemented for complex wi')
end

arg = ob.arg;
arg.mask = ob.imask_array;
[arg.fftkern arg.reuse] = Gdsft_gram_init(arg, wi, reuse, ob.idim);

if any(~arg.mask(:))
	omask = arg.mask;
else
	omask = []; % trick: omask not needed in full case
end
T = fatrix2('arg', arg', ...
	'idim', ob.idim, 'imask', arg.mask, ...
	'odim', ob.idim, 'omask', omask, ...
	'forw', forw, 'back', back);


% Gdsft_gram_init()
% construct kernel of circulant matrix into which T is embedded
% and take its DFT to prepare for multiplication.
% out
%	fftkern [[2Nd]]		FFT of kernel of circulant matrix
%
function [fftkern, reuse] = Gdsft_gram_init(arg, wi, reuse, Nd)

switch numel(Nd)
case 1
	[fftkern reuse] = Gdsft_gram_init1(arg, wi, reuse, Nd, arg.show);
case 2
	[fftkern reuse] = Gdsft_gram_init2(arg, wi, reuse, Nd, arg.show);
case 3
	[fftkern reuse] = Gdsft_gram_init3(arg, wi, reuse, Nd, arg.show);
otherwise
	fail('dim %d not done', numel(Nd))
end


% Gdsft_gram_init1()
% 1d filter for circulant matrix
% note: only toeplitz kernel values from -(N-1) to (N-1) are relevant
% so the value at +/- N does not matter so we set it to zero.
function [fftkern, reuse] = Gdsft_gram_init1(arg, wi, reuse, N1, show)

if isempty(reuse)
	reuse.G1 = Gdsft(arg.om_t', arg.Nd); % with n_shift = 0
end

block1 = reuse.G1' * real(wi); % kludge

% kernel of Toeplitz matrix from -N to N-1 but with fftshift
% This is inherently Hermitian symmetric except for the middle [0] value.
% Use 0 for the irrelevant value at -N.
err1 = abs(imag(block1(1))) / abs(block1(1));
tol = 0;
if err1 > tol
	printm('removing imaginary h[0] part of relative size %g', err1)
	block1(1) = real(block1(1));
end
kern = [block1; 0; flipud(conj(block1(2:N1)))]; % [2*N1]
if show
	kern = dsft_gram_hermitify(kern, show); % need not due to block1(1) fix
end
fftkern = fft(kern); % [2*N1]


% Gdsft_gram_init2()
function [fftkern, reuse] = Gdsft_gram_init2(arg, wi, reuse, Nd, show)

if isempty(reuse)
	[reuse.G1 reuse.G2] = Gdsft_gram_setup2(arg);
end

N1 = Nd(1);
N2 = Nd(2);

block1 = reshape(reuse.G1' * real(wi), Nd); % kludge
block2 = reshape(reuse.G2' * real(wi), Nd);

z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
	[block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
	zeros(1,2*N2);
	[flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
]; % [(2Nd)]
kern = dsft_gram_hermitify(kern, show); % force Hermitian symmetry
fftkern = fftn_fast(kern);


% Gdsft_gram_init3()
function [fftkern reuse] = Gdsft_gram_init3(arg, wi, reuse, Nd, show)
if isempty(reuse)
	[reuse.G1 reuse.G2 reuse.G3 reuse.G4 ] = build_G1_G2_G3_G4(arg);
end
N1 = Nd(1);
N2 = Nd(2);
N3 = Nd(3);

filtblk1 = reshape(reuse.G1' * real(wi), Nd);
filtblk2 = reshape(reuse.G2' * real(wi), Nd);
filtblk3 = reshape(reuse.G3' * real(wi), Nd);
filtblk4 = reshape(reuse.G4' * real(wi), Nd);

tblk1 = filtblk1;
tblk2 = filtblk2(2:N1,:,:);	% remove the duplicated part with filtblk1
tblk3 = filtblk3(2:N1,2:N2,:);	% remove the duplicated part with block 1, 2, 4
tblk4 = filtblk4(:,2:N2,:);	% remove the duplicated part with block 1

z1 = zeros(N1,1,N3);
z2 = zeros(N1-1,1,N3);

% top half of the 3D filter
kern_top = [
	tblk1 z1 flipdim(tblk4,2);	% upper block
	zeros(1,2*N2,N3);		% zero padding in the middle
	flipdim(tblk2,1) z2 flipdim(flipdim(tblk3,1),2) % lower block
];

% construct the bottom half now
bblk1 = flipdim(conj(filtblk3(:,:,2:N3)),3);
bblk2 = flipdim(conj(filtblk4(:,:,2:N3)),3);
bblk2 = bblk2(2:N1,:,:);
bblk3 = flipdim(conj(filtblk1(:,:,2:N3)),3);
bblk3 = bblk3(2:N1,2:N2,:);
bblk4 = flipdim(conj(filtblk2(:,:,2:N3)),3);
bblk4 = bblk4(:,2:N2,:);

z4 = zeros(N1,1,N3-1);
z5 = zeros(N1-1,1,N3-1);
kern_bottom = [
	bblk1 z4 flipdim(bblk4,2);
	zeros(1,2*N2,N3-1);
	flipdim(bblk2,1) z5 flipdim(flipdim(bblk3,1),2)
];

kern = cat(3, kern_top, zeros(2*N1,2*N2,1));
kern = cat(3, kern, kern_bottom);

kern = dsft_gram_hermitify(kern, show); % force Hermitian symmetry
fftkern = fftn_fast(kern);


% Gdsft_gram_setup2()
% modified versions of A for 2D
% this routine is quite tricky.  nothing new is really rebuilt here,
% but some structure values are changed.
function [G1, G2] = Gdsft_gram_setup2(arg)
om = arg.om_t';
Nd = arg.Nd;
G1 = Gdsft(om, Nd);
om(:,1) = -om(:,1);
G2 = Gdsft(om, Nd);


% build_G1_G2_G3_G4()
% Created by Daehyun Yoon for 3D case
% build more Gnufft objects, with negative om.
function [G1, G2, G3, G4] = build_G1_G2_G3_G4(arg)
om = arg.om_t';
Nd = arg.Nd;
G1 = Gdsft(om, Nd);
om(:,1) = -om(:,1);
G2 = Gdsft(om, Nd);
om(:,2) = -om(:,2);
G3 = Gdsft(om, Nd);
om(:,1) = -om(:,1);
G4 = Gdsft(om, Nd);


% Gdsft_gram_mult()
% multiply an image x by a toeplitz matrix
% by embedding it into a circulant matrix of twice the size and using FFT.
% in
%	x	[(Nd)]
%	fftkern	[[2Nd]]
%	mask	[(Nd)]
% out
%	y	[(Nd)]
function y = Gdsft_gram_mult(x, fftkern, mask)

N2 = size(fftkern);
if numel(N2) == 2 && N2(2) == 1
	N2 = N2(1);
end
Nd = N2 / 2;

LL = size(x,1+numel(N2));
y = zeros([Nd LL]);

switch numel(N2)
case 1
	tmp = fft(x, N2); % [N2 L]
	tmp = tmp .* repmat(fftkern, [1 ncol(tmp)]);
	y = ifft(tmp);
	y = y(1:Nd,:); % [Nd L]
	y = y .* repmat(mask, [1 ncol(y)]); % fatrix2 requires mask

case 2
	for ll=1:LL
		tmp = ifftn_fast(fftkern .* fftn_fast(x(:,:,ll), N2));
		tmp = tmp(1:Nd(1), 1:Nd(2));
		tmp = tmp .* mask; % fatrix2 requires mask
		y(:,:,ll) = tmp;
	end

case 3
	for ll=1:LL
		tmp = ifftn_fast(fftkern .* fftn_fast(x(:,:,:,ll), N2));
		tmp = tmp(1:Nd(1), 1:Nd(2), 1:Nd(3));
		tmp = tmp .* mask; % fatrix2 requires mask
		y(:,:,:,ll) = tmp;
	end

otherwise
	fail('dim %d not done', numel(N2))
end



% Gdsft_gram_test1()
% test 1d version
function Gdsft_gram_test1
N = 16;
M = 21;
omega = 2*pi*rand(M,1);
wi = [1:size(omega,1)]';
mask = true(N,1); mask(end-3) = false;
A = Gdsft(omega, N, 'mask', mask, 'n_shift', N/2, 'show', 1);
T = build_gram(A, wi);
test_adjoint(T, 'complex', 1, 'tolre', 1e-7);
fatrix2_tests(A, 'complex', 1, 'tol_gram', 2e-7);


% Gdsft_gram_test2()
function Gdsft_gram_test2
N = [16 14];
fov = N;
ktype = 'spiral0';
[kspace omega wi] = mri_trajectory(ktype, {}, N, fov, {'voronoi'});
ig = image_geom('nx', N(1), 'ny', N(2), 'dx', 1);
mask = ellipse_im(ig, [0 0 14 15 0 1], 'oversample', 3) > 0;
A = Gdsft(omega, N, 'n_shift', N/2, 'show', 1);
fatrix2_tests(A, 'complex', 1, 'tol_gram', 2e-7);

T = build_gram(A, wi);
fatrix2_tests(T, 'complex', 1, 'tol_gram', 1e-7);
test_adjoint(T, 'complex', 1, 'tolre', 1e-7); % matches!

if 1
	x = ellipse_im(ig, 'shepplogan-emis');
	x1 = A' * (wi .* (A * x));
	x1 = embed(x1, mask);
	x2 = T * x;
	equivs(x1, x2, 'thresh', 2e-7)
end

tic
for ii=1:50, b1 = embed(A' * (wi .* (A * x(mask))), mask); end
t1 = toc;
tic
for ii=1:50, b2 = embed(T * x(mask), mask); end
t2 = toc;
equivs(b1, b2)
printm('time: A''Ax = %.3f, Tx=%.3f', t1, t2)


% Gdsft_gram_test3()
function Gdsft_gram_test3
N = [8 4 2];
rng(0)
omega = rand(21,3) * 2 * pi;
A = Gdsft(omega, N, 'n_shift', N/2, 'show', 1);
fatrix2_tests(A, 'complex', 1, 'tol_gram', 2e-7);
test_adjoint(A, 'complex', 1, 'big', 1);

wi = [1:size(omega,1)]';
T = build_gram(A, wi);
fatrix2_tests(T, 'complex', 1, 'tol_gram', 1e-7);
test_adjoint(T, 'complex', 1, 'tolre', 1e-7);

if 1
	x = cumsum(ones(N), 2);
	x1 = A' * (wi .* (A * x));
	x1 = iembed(A, x1);
	x2 = T * x;
	equivs(x1, x2)
end


% Gdsft_gram_test()
function Gdsft_gram_test
Gdsft_gram_test1
Gdsft_gram_test2
Gdsft_gram_test3
