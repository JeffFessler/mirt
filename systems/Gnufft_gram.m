  function [T, reuse] = Gnufft_gram(ob, W, reuse)
%|function [T, reuse] = Gnufft_gram(A, W, reuse)
%|
%| construct Toeplitz gram matrix object T = A'WA for a Gnufft object
%| in
%|	A	[M np]		Gnufft object (Fatrix or fatrix2)
%|	W	[M M]		W = diag_sp(wi) (often simply "1" or [])
%|				W = Gdiag() for fatrix2
%|	reuse	struct		stuff from the last call that can be reused
%|
%| out
%|	T	[np np]		fatrix2 or Fatrix object
%|	reuse	struct		stuff that can be reused on future calls
%|
%| Copyright 2004-6-29, Jeff Fessler & Hugo Shi, University of Michigan

if nargin == 1 && streq(ob, 'test'), Gnufft_gram_test, return, end
if nargin == 1 && streq(ob, 'test1'), Gnufft_gram_test1, return, end
if nargin == 1 && streq(ob, 'test2'), Gnufft_gram_test2, return, end
if nargin < 2 || nargin > 3, help(mfilename), error(mfilename), end

if ~isvar('reuse'), reuse = []; end

if ~isa(ob, 'fatrix2') && ~isa(ob, 'Fatrix')
	fail('A wrong class: %s', class(ob))
end

[nd np] = size(ob);
forw = @(arg, x) Gnufft_gram_mult(x, ...
	arg.fftkern, arg.mask, np, streq(class(ob), 'Fatrix'));

% extract wi
if isnumeric(W) % includes case where W is empty
	wi = W;
elseif isa(W, 'Fatrix') && streq(W.caller, 'diag_sp')
	wi = W.arg.diag;
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
	fail('not implemented for complex wi (see old version in arch/)')
end

arg = ob.arg;

switch class(ob)
case 'Fatrix'
	Nd = size(arg.mask);
	if numel(Nd) == 2 && Nd(2) == 1
		Nd = Nd(1); % trick for 1d
	end
	[arg.fftkern arg.reuse] = ...
		Gnufft_gram_init(arg, wi, reuse, np, Nd);
	T = Fatrix([1 1]*np, arg, 'forw', forw, 'back', back);

case 'fatrix2'
	arg.mask = ob.imask_array;
	[arg.fftkern arg.reuse] = ...
		Gnufft_gram_init(arg, wi, reuse, np, ob.idim);

	if any(~arg.mask(:))
		omask = arg.mask;
	else
		omask = []; % trick: omask not needed in full case
	end
	T = fatrix2('arg', arg', ...
		'idim', ob.idim, 'imask', arg.mask, ...
		'odim', ob.idim, 'omask', omask, ...
		'forw', forw, 'back', back);

otherwise
	fail('class "%s" unknown', class(ob))
end


% Gnufft_gram_init()
% construct kernel of circulant matrix into which T is embedded
% and take its DFT to prepare for multiplication.
% out:
% fftkern [[2Nd]]	FFT of kernel of circulant matrix
%
function [fftkern, reuse] = Gnufft_gram_init(arg, wi, reuse, np, Nd)

if isfield(arg, 'show')
	show = arg.show;
else
	show = false;
end
switch numel(Nd)
case 1
	[fftkern reuse] = Gnufft_gram_init1(arg, wi, reuse, np, Nd, show);
case 2
	[fftkern reuse] = Gnufft_gram_init2(arg, wi, reuse, np, Nd, show);
case 3
	[fftkern reuse] = Gnufft_gram_init3(arg, wi, reuse, np, Nd, show);
otherwise
	fail('dim %d not done', numel(Nd))
end


% Gnufft_gram_init1()
% 1d filter for circulant matrix
% note: only toeplitz kernel values from -(N-1) to (N-1) are relevant
% so the value at +/- N does not matter so we set it to zero.
function [fftkern, reuse] = Gnufft_gram_init1(arg, wi, reuse, np, N1, show)

if isempty(reuse)
	reuse.G1 = Gnufft_gram_setup1(arg);
end

block1 = reuse.G1' * real(wi); % kludge

% kernel of Toeplitz matrix from -N to N-1 but with fftshift
% this is inherently Hermitian symmetric except for the middle [0] value
% use 0 for the irrelevant value at -N
err1 = abs(imag(block1(1))) / abs(block1(1));
tol = 0;
if err1 > tol
	printm('removing imaginary h[0] part of relative size %g', err1)
	block1(1) = real(block1(1));
end
kern = [block1; 0; flipud(conj(block1(2:N1)))]; % [2*N1]
kern = dsft_gram_hermitify(kern, show); % force Hermitian
fftkern = fft(kern); % [2*N1]


% Gnufft_gram_init2()
function [fftkern, reuse] = Gnufft_gram_init2(arg, wi, reuse, np, Nd, show)
if 0 % old
	g1 = build_Gmod1(arg);
	g2 = build_Gmod2(arg);
	s1 = g1.arg.st;
	s2 = g2.arg.st;
%	max_percent_diff(s1.phase_shift, s2.phase_shift)
	printm('scaling factors: mpd=%g%%', max_percent_diff(s1.sn, s2.sn))
%	max_percent_diff(s1.alpha{1}, s2.alpha{1})
	printm('table1: mpd=%g%%', max_percent_diff(s1.h{1}, s2.h{1}))
end

if 0 % old
	kern = fftshift(reshape(g2' * wi, g2.arg.st.Nd));
	fftkern = fftn_fast(kern);
return
	t = kern;
	%t = t ./ Gmod.st.sn;
	t(end/2+1,:) = 0;
	t(:,end/2+1) = 0;
end

if isempty(reuse)
	[reuse.G1 reuse.G2] = Gnufft_gram_setup2(arg);
end

N1 = Nd(1);
N2 = Nd(2);

block1 = reshape(reuse.G1' * real(wi), Nd); % kludge
block2 = reshape(reuse.G2' * real(wi), Nd);

if 0 % hugo way
	kern = zeros(2*Nd);
	kern(1:N1,1:N2) = block1;
	kern(N1+2:2*N1, N2+2:2*N2) = ...
		fliplr(flipud(conj(block1(2:N1, 2:N2))));
	kern(N1+2:2*N1, 1:N2) = flipud(block2(2:N1, :));
	kern(1:N1, N2+2:2*N2) = ...
		conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]));
end

z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
	[block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
	zeros(1,2*N2);
	[flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
]; % [(2Nd)]
kern = dsft_gram_hermitify(kern, show); % force Hermitian symmetry
fftkern = fftn_fast(kern);


% Gnufft_gram_init3()
function [fftkern reuse] = Gnufft_gram_init3(arg, wi, reuse, np, Nd, show)
if isempty(reuse)
	[reuse.G1 reuse.G2 reuse.G3 reuse.G4 ] = build_G1_G2_G3_G4(arg);
end
N1 = Nd(1);
N2 = Nd(2);
N3 = Nd(3);

Filtblk1 = reshape(reuse.G1' * real(wi), Nd);
Filtblk2 = reshape(reuse.G2' * real(wi), Nd);
Filtblk3 = reshape(reuse.G3' * real(wi), Nd);
Filtblk4 = reshape(reuse.G4' * real(wi), Nd);

tblk1 = Filtblk1;
tblk2 = Filtblk2(2:N1,:,:);	%remove the duplicated part with Filtblk1
tblk3 = Filtblk3(2:N1,2:N2,:);	%remove the duplicated part with block 1, 2, 4
tblk4 = Filtblk4(:,2:N2,:);	%remove the duplicated part with block 1

z1 = zeros(N1,1,N3);
z2 = zeros(N1-1,1,N3);

% top half of the 3D filter
kern_top = [
	tblk1 z1 flipdim(tblk4,2);	% Upper block
	zeros(1,2*N2,N3);		% Zero padding in the middle
	flipdim(tblk2,1) z2 flipdim(flipdim(tblk3,1),2) % lower block
];

% construct the bottom half now
bblk1 = flipdim(conj(Filtblk3(:,:,2:N3)),3);
bblk2 = flipdim(conj(Filtblk4(:,:,2:N3)),3);
bblk2 = bblk2(2:N1,:,:);
bblk3 = flipdim(conj(Filtblk1(:,:,2:N3)),3);
bblk3 = bblk3(2:N1,2:N2,:);
bblk4 = flipdim(conj(Filtblk2(:,:,2:N3)),3);
bblk4 = bblk4(:,2:N2,:);

z4 = zeros(N1,1,N3-1);
z5 = zeros(N1-1,1,N3-1);
kern_bottom = [
	bblk1 z4 flipdim(bblk4,2);
	zeros(1,2*N2,N3-1);
	flipdim(bblk2,1) z5 flipdim(flipdim(bblk3,1),2)
];

kern = cat(3,kern_top,zeros(2*N1,2*N2,1));
kern = cat(3,kern,kern_bottom);

kern = dsft_gram_hermitify(kern, show); % force Hermitian symmetry
fftkern = fftn_fast(kern);


% Gnufft_gram_setup1()
% modified version of G for 1D case
% this routine is quite tricky.  nothing new is really rebuilt here,
% but some structure values are changed.
function G1 = Gnufft_gram_setup1(arg)
st = arg.st;
mask = true(st.Nd, 1); % trick: full mask!

if isfield(st, 'phase_shift') % table-based nufft
	st = rmfield(st, 'phase_shift'); % trick: eliminate phase shift
	st.n_shift = 0 * st.n_shift;
	G1 = Gnufft(mask, st);
else
%	warning 'non-table nufft may not work'
	G1 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0]});
end


% Gnufft_gram_setup2()
% modified versions of G for 2D
% this routine is quite tricky.  nothing new is really rebuilt here,
% but some structure values are changed.
function [G1, G2] = Gnufft_gram_setup2(arg)
st = arg.st;
mask = true(st.Nd); % trick: full mask!

if isfield(st, 'phase_shift') % table-based nufft
	st = rmfield(st, 'phase_shift'); % trick: eliminate phase shift
	st.n_shift = 0 * st.n_shift;

	% build two Gnufft objects, one with negative om1.
	G1 = Gnufft(mask, st);
	st.om(:,1) = -st.om(:,1);
	G2 = Gnufft(mask, st);

else
%	warning 'non-table nufft may not work'
	G1 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0]});
	st.om(:,1) = -st.om(:,1);
	G2 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0]});
end


% build_G1_G2_G3_G4()
% Created by Daehyun Yoon for 3D case
function [G1, G2, G3, G4] = build_G1_G2_G3_G4(arg)
st = arg.st;
mask = true(st.Nd); % trick: full mask!

if isfield(st, 'phase_shift') % table-based nufft
	st = rmfield(st, 'phase_shift'); % trick: eliminate phase shift
	st.n_shift = 0 * st.n_shift;

	% now build more Gnufft objects, with negative om.
	G1 = Gnufft(mask, st);
	st.om(:,1) = -st.om(:,1);
	G2 = Gnufft(mask, st);
	st.om(:,2) = -st.om(:,2);
	G3 = Gnufft(mask, st);
	st.om(:,1) = -st.om(:,1);
	G4 = Gnufft(mask,st);

else
%	warning 'non-table nufft may not work'
	G1 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0 0]});
	st.om(:,1) = -st.om(:,1);
	G2 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0 0]});
	st.om(:,2) = -st.om(:,2);
	G3 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0 0]});
	st.om(:,1) = -st.om(:,1);
	G4 = Gnufft(mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0 0]});
end


% Gnufft_gram_mult()
% multiply an image x by a toeplitz matrix
% by embedding it into a circulant matrix of twice the size and using FFT.
% in
%	x	[*Nd 1] or [np 1]
%	fftkern	[[2Nd]]
%	mask	[[Nd]]
% out
%	y	[*Nd 1] or [np 1]
function y = Gnufft_gram_mult(x, fftkern, mask, np, do_Fatrix)

N2 = size(fftkern);
if numel(N2) == 2 && N2(2) == 1
	N2 = N2(1);
end
Nd = N2 / 2;

if size(x,1) == np
	x = embed(x, mask); % [np (R)] to [(Nd) (R)]
end

LL = size(x,1+numel(N2));
y = zeros([Nd LL]);

switch numel(N2)
case 1
	tmp = fft(x, N2); % [N2 L]
	tmp = tmp .* repmat(fftkern, [1 ncol(tmp)]);
	y = ifft(tmp);
	y = y(1:Nd,:); % [Nd L]
	if ~do_Fatrix % fatrix2 requires mask
		y = y .* repmat(mask, [1 ncol(y)]);
	end

case 2
	for ll=1:LL
		tmp = ifftn_fast(fftkern .* fftn_fast(x(:,:,ll), N2));
		tmp = tmp(1:Nd(1), 1:Nd(2));
		if ~do_Fatrix % fatrix2 requires mask
			tmp = tmp .* mask;
		end
		y(:,:,ll) = tmp;
	end

case 3
	for ll=1:LL
		tmp = ifftn_fast(fftkern .* fftn_fast(x(:,:,:,ll), N2));
		tmp = tmp(1:Nd(1), 1:Nd(2), 1:Nd(3));
		if ~do_Fatrix % fatrix2 requires mask
			tmp = tmp .* mask;
		end
		y(:,:,:,ll) = tmp;
	end

otherwise
	fail('dim %d not done', numel(N2))
end

if do_Fatrix % for Fatrix, must always make output column
	y = reshapee(y, [], LL);
	y = y(mask(:),:);
else % for fatrix2
%	y = y .* repmat(mask, [ones(1,ncol(y)-1) LL]); % apply mask
end


% Gnufft_gram_test1()
% test 1d version
function Gnufft_gram_test1
N = 16;
J = 6;
K = 2*N;
rng(0);
M = 51;
omega = 2*pi*rand(M,1);
nufft_args = {N, J, K, N/2, 'table', 2^11, 'minmax:kb'};
wi = [1:size(omega,1)]';
mask = true(N,1); mask(end-3) = false;
A = Gnufft('fatrix2', mask, {omega, nufft_args{:}});
T = build_gram(A, wi);
test_adjoint(T, 'complex', 1, 'tolre', 1e-7);
fatrix2_tests(A, 'complex', 1, 'tol_gram', 2e-5); % needed due to NUFFT approx


% Gnufft_gram_test2()
% test this object and compare its speed to Gnufft approach
function Gnufft_gram_test2
N = [16 14];
fov = N;
J = [6 5+0];
K = 2*N;
ktype = 'spiral0';
[kspace omega wi] = mri_trajectory(ktype, {}, N, fov, {'voronoi'});
nufft_args = {N, J, K, N/2, 'table', 2^12, 'minmax:kb'};
ig = image_geom('nx', N(1), 'ny', N(2), 'dx', 1);
mask = ellipse_im(ig, [0 0 14 15 0 1], 'oversample', 3) > 0;
x = ellipse_im(ig, 'shepplogan-emis');

if 1 % compare to exact Gdsft
	A = Gnufft(mask, {omega, nufft_args{:}});
	G = Gdsft(omega, N, 'n_shift', N/2);
	A = A(:,:);
	G = G(:,:);
	im plc 2 3
	im(1, real(A)')
	im(2, real(G)')
	im(3, real(A - G)')
	im(4, imag(A)')
	im(5, imag(G)')
	im(6, imag(A - G)')
	equivs(A, G, 'thresh', 8e-4) % surprisingly large threshold needed?
prompt
end

classes = {'fatrix2', 'Fatrix'};
for ic=1:numel(classes)
	pr ic
	A = Gnufft(classes{ic}, mask, {omega, nufft_args{:}});
	fatrix2_tests(A, 'complex', 1, 'tol_gram', 3e-5);

	T = build_gram(A, diag_sp(wi));
	fatrix2_tests(T, 'complex', 1);

	tic
	for ii=1:50, b1 = embed(A' * (wi .* (A * x(mask))), mask); end
	t1 = toc;
	tic
	for ii=1:50, b2 = embed(T * x(mask), mask); end
	t2 = toc;
	printm('time: A''Ax = %.3f, Tx=%.3f', t1, t2)

	equivs(b1, b2, 'thresh', 7e-5) % bigger threshold due to NUFFT approx
	if 0
		d = max_percent_diff(b1, b2);
		printm('max percent diff between Tx and A''Ax = %g%%', d)
		im plc 2 1
		im(1, abs(stackup(b1, b2)), 'A''A and T'), cbar
		im(2, abs(b1-b2), 'diff'), cbar
	end
end


% Gnufft_gram_test()
function Gnufft_gram_test()
Gnufft_gram_test1
Gnufft_gram_test2
