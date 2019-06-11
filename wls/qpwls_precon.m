 function M = qpwls_precon(type, sys, C, mask, varargin)
%function M = qpwls_precon(type, sys, C, mask, varargin)
%|
%| build a preconditioner for QPWLS problems:
%| approximations to inv([A'WA + C'C])
%|
%| usage:
%| M = qpwls_precon('circ0', {T}, C, mask);
%| M = qpwls_precon('circ0', {A, W}, C, mask);
%|
%| in
%|	type	string		'circ0' : circulant based on center of FOV
%|				'dcd0' : diagonal * circulant * diagonal
%|	sys	cell		{T = A'WA} or {A, W}
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix
%|	C	[nc np]		penalty 'derivatives' (R = \Half C'*C)
%|	mask	[nx ny]		which pixels are updated
%|
%| options
%|	kappa	[nx ny]		needed for dcd0
%|	apod	[nx ny]		optional apodizer for A'WAe
%|	'class'	''		'fatrix2' or 'Fatrix' (default)
%|
%| out
%|	M	[np np]		fatrix2 (default) or Fatrix object
%|
%| The 'dcd0' preconditioner is based on Fessler&Booth IEEE T-IP May 1999
%|
%| Copyright 2004-6-29, Jeff Fessler, University of Michigan

if nargin == 1 && streq(type, 'test'), qpwls_precon_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

arg.chat = false;
arg.kappa = [];
arg.apod = [];
arg.class = 'fatrix2';
arg = vararg_pair(arg, varargin);

arg.mask = mask;
idim = size(mask);

switch type
case 'circ0' % circulant preconditioner
	arg.Mdft = qpwls_precon_circ_Mdft(sys, C, mask, arg.apod, arg.chat);
	arg.st = strum(arg, {'plot', @qpwls_precon_plot, '()'});
	switch arg.class
	case 'Fatrix'
		dim = [1 1] * numel(mask);
		M = Fatrix(dim, arg, 'forw', @qpwls_precon_circ_mult_Fatrix);
	case 'fatrix2'
		M = fatrix2('idim', idim, 'odim', idim, 'arg', arg, ...
			'imask', mask, 'omask', mask, ...
			'forw', @qpwls_precon_circ_mult);
	otherwise
		fail 'class'
	end

case 'dcd0' % diagonal / circulant / diagonal
	arg = qpwls_precon_dcd0_init(sys, C, arg);
	switch arg.class
	case 'Fatrix'
		dim = [1 1] * numel(mask);
		M = Fatrix(dim, arg, 'forw', @qpwls_precon_dcd0_mult_Fatrix);
	case 'fatrix2'
		M = fatrix2('idim', idim, 'odim', idim, 'arg', arg, ...
			'imask', mask, 'omask', mask, ...
			'forw', @qpwls_precon_dcd0_mult);
	otherwise
		fail 'class'
	end

case 'diag' % diagonal
	fail 'todo: ask jeff'

otherwise
	fail('unknown preconditioner "%s"', type)
end


% qpwls_precon_plot()
% show spectrum and psf
function psf = qpwls_precon_plot(st)

[nx ny] = size(st.mask);
tmp = zeros(nx, ny);
tmp(nx/2+1,ny/2+1) = 1;
ix = [-nx/2:nx/2-1];
iy = [-ny/2:ny/2-1];

im plc 1 2
im(1, ix, iy, fftshift(st.Mdft)), cbar
xtick([-nx/2 0 nx/2-1])
ytick([-ny/2 0 ny/2-1])

psf = qpwls_precon_circ_mult_Fatrix(st, tmp(st.mask));
psf = embed(psf, st.mask);
%psf = reale(fftshift(ifftn(st.Mdft))); % same
im(2, ix, iy, psf), cbar
xtick([-nx/2 0 nx/2-1])
ytick([-ny/2 0 ny/2-1])


% qpwls_precon_dcd0_init()
% initialize dcd preconditioner based on center pixel
function arg = qpwls_precon_dcd0_init(sys, R, arg)
if length(sys) ~= 2, error 'need cell(2)', end
switch arg.class
case 'Fatrix'
	arg.diag = 1 ./ arg.kappa(arg.mask(:));
case 'fatrix2'
	arg.diag = div0(1, arg.kappa);
end

ej = qpwls_precon_e0(arg.mask);
ctc = R.cgrad(R, 1e-2 * ej) / 1e-2; % trick
ctc = embed(ctc, arg.mask);
ctc = ctc / arg.kappa(end/2+1,end/2+1)^2; % trick

sys = {sys{1}, 1}; % trick
arg.Mdft = qpwls_precon_circ_Mdft(sys, ctc, arg.mask, [], arg.chat);


% qpwls_precon_dcd0_mult_Fatrix()
% multiply using diag * circ * diag
function y = qpwls_precon_dcd0_mult_Fatrix(arg, x)

x = x .* arg.diag;
y = ifftn_fast(arg.Mdft .* fftn_fast(embed(x, arg.mask)));
y = y(arg.mask(:));
y = y .* arg.diag;
if isreal(x)
	y = reale(y, 'warn');
end


% qpwls_precon_dcd0_mult()
% multiply using diag * circ * diag
function y = qpwls_precon_dcd0_mult(arg, x)

x = x .* arg.diag;
y = ifftn_fast(arg.Mdft .* fftn_fast(x));
y = y .* arg.diag;
if isreal(x)
	y = reale(y, 'warn');
end


% qpwls_precon_circ_Mdft()
% setup circulant preconditioner based on center pixel
function Mdft = qpwls_precon_circ_Mdft(sys, C, mask, apod, chat)

ej = qpwls_precon_e0(mask);

% T * x or A'WA*x
if ~iscell(sys), error 'sys must be cell', end
if length(sys) == 2
	A = sys{1};
	W = sys{2};
	awa = A' * (W * (A * ej));
elseif length(sys) == 1
	awa = sys{1} * ej; % T * ej
else
	error 'unknown cell'
end

awa = embed(awa, mask);
if ~isempty(apod)
%	im(awa, 'awa'), cbar, prompt
	awa = awa .* apod;
%	im(awa, 'awa apodized'), cbar, prompt
end

if isnumeric(C) && isequal(size(C), size(awa)) % trick
	ccc = C;

elseif isstruct(C) % 'R'
	R = C;
	ccc = R.cgrad(R, 1e-2 * ej) / 1e-2; % trick
	ccc = embed(ccc, mask);

else
	ccc = C' * (C * ej);
	ccc = embed(ccc, mask);
end
%im(ccc, 'ccc'), cbar, prompt

f.awa = fftn_fast(ifftshift(awa));
f.ccc = fftn_fast(ifftshift(ccc));
f.ccc = reale(f.ccc, 'warn'); % these should be nearly real
if any(f.ccc(:) < - 1e-6 * max(f.ccc(:)))
	printm('ccc min = %g max = %g', min(f.ccc(:)), max(f.ccc(:)))
	clf, im(ccc), keyboard
	error 'bug: circulant penalty is not nonnegative definite!?'
end
f.ccc = max(f.ccc, 0);
f.awa = reale(f.awa, 'warn');
if min(f.awa(:)) < 0
	printm('setting %g%% to zero', ...
		min(f.awa(:)) / max(f.awa(:)) * 100)
	f.awa = max(f.awa, 0);
end
%minmax(f.awa), minmax(f.ccc)
f.h = f.awa + f.ccc;	% approximate hessian in freq. domain
if min(f.h(:)) <= 0
	error 'circulant preconditioner has zero!?  do you regularize?'
end
if chat
	printm('approximate condition number: %g', ...
		max(f.h(:)) / min(f.h(:)))
end

Mdft = 1 ./ f.h;



% qpwls_precon_circ_mult_Fatrix()
% multiply using fft
function y = qpwls_precon_circ_mult_Fatrix(arg, x)

y = ifftn_fast(arg.Mdft .* fftn_fast(embed(x, arg.mask)));
y = y(arg.mask(:));
if isreal(x)
	y = reale(y, 'warn');
end


% qpwls_precon_circ_mult()
% multiply using fft
function y = qpwls_precon_circ_mult(arg, x)

y = ifftn_fast(arg.Mdft .* fftn_fast(x));
%y = y(arg.mask(:));
if isreal(x)
	y = reale(y, 'warn');
end


% qpwls_precon_e0()
% unit impulse at "center" voxel
function e0 = qpwls_precon_e0(mask)

e0 = zeros(size(mask));
switch ndims(mask)
case 2
	if size(mask,2) == 1 % 1d
		e0(end/2+1) = 1;
	else % 2d
		e0(end/2+1,end/2+1) = 1;
	end
case 3
	e0(end/2+1,end/2+1,end/2+1) = 1;
otherwise
	fail 'not done'
end
e0 = e0(mask(:));


% qpwls_precon_test_circ0()
function qpwls_precon_test_circ0(mask, A, reg, type_diff)

im(2, mask)

N = size(mask);
C = Cdiffs(N, 'offsets', [1 N(1)], 'mask', mask, 'type_diff', type_diff);
C = sqrt(reg) * C;

M = qpwls_precon('circ0', {A, 1}, C, mask); % preconditioner

Af = full(A);
Cf = full(C);
H = Af'*Af + Cf'*Cf;

im(3, H)
titlef('cond(H) = %g', cond(H))

Mf = full(M);
% Mf = M * eye(sum(mask(:)));
im(4, Mf);

MH = Mf * H;
im(5, MH)
titlef('cond(MH) = %g', cond(MH))


% qpwls_precon_test
function qpwls_precon_test

im plc 2 3

if 0 % test deconvolution problems
	N = [16 14];
	mask = true(N);
	psf = ones(5)/5^2;
	% periodic case where circ0 is perfect
	A = Gblur(mask, 'psf', psf, 'type', 'conv,per'); % easy
	im(1, psf)
	qpwls_precon_test_circ0(mask, A, 2^1, 'circshift'); % periodic
	prompt

	% non-periodic case where circ0 is approximate
	mask(1) = false; % stress
	A = Gblur(mask, 'psf', psf, 'type', 'conv,same');
	qpwls_precon_test_circ0(mask, A, 2^1, '');
return
end

if 1 % NUFFT case
	N = [32 30];
	N = [16 14]; % todo
	J = [6 6];
	K = 2*N;
	% spiral trajectory with reasonable sampling along k-space axes.
	omega = linspace(0, (max(N)-1)*2*pi/2, max(N)^2+1)';
	omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);

	if im
		im subplot 1, plot(omega(:,1), omega(:,2), '.')
		axis_pipi, axis equal
		titlef('%d k-space samples', size(omega,1))
	end

	mask = true(N);
	mask(1) = false; % stress

	A = Gnufft(mask, {omega, N, J, K});
	qpwls_precon_test_circ0(mask, A, 2^2, 'circshift');
end
