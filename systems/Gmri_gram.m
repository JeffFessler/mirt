  function [T, reuse] = Gmri_gram(ob, W, reuse)
%|function [T, reuse] = Gmri_gram(ob, W, reuse)
%|
%| build Toeplitz-like gram-matrix object for G'WG gram matrix
%| to be called by build_gram() indirectly rather than directly by user!
%|
%| Jeff Fessler

if nargin < 1, help(mfilename), error(mfilename), end
if streq(ob, 'test'), Gmri_gram_test, return, end

if isempty(W)
	wi = 1; % default unweighted

% 1D column vector or scalar:
elseif isnumeric(W) && ndims(W) == 2 && size(W,2) == 1
	wi = W;
elseif isa(W, 'Fatrix') && streq(W.caller, 'diag_sp')
	wi = W.arg.diag;
elseif isa(W, 'fatrix2') && streq(W.caller, 'Gdiag')
	wi = W.arg.diag;
else
	fail('Gmri_gram requires W to be Gdiag or diag_sp or wi array')
end

if ~isreal(wi), fail('only real wi supported; otherwise T not Hermitian'), end

arg2.wi = wi;
arg2.Gmri = ob;
arg2.new_zmap = @Gmri_gram_new_zmap;

T = Gmri_gram_work(arg2, reuse);


% Gmri_gram_work()
function T = Gmri_gram_work(arg2, reuse)

arg1 = arg2.Gmri.arg;
wi = arg2.wi;

if isempty(arg1.zmap)
	T = build_gram(arg1.Gnufft, wi .* abs(arg1.basis.transform).^2);
else
	LL = ncol(arg1.aB);
	arg2.T = cell(LL,1);
	reuse = [];
	for ll=1:LL
		wl = arg1.aB(:,ll) .* wi .* abs(arg1.basis.transform).^2;
		[arg2.T{ll} reuse] = build_gram(arg1.Gnufft, wl, reuse);
	end

	arg2.dim = arg1.dim([2 2]); % [np np]

	switch class(arg2.Gmri)
	case 'Fatrix'
		T = Fatrix(arg2.dim, arg2, ...
			'caller', [mfilename '.Gmri_gram'], ...
			'forw', @Gmri_zmap_forw_Fatrix, ...
			'back', @Gmri_zmap_forw_Fatrix); % trick: Hermitian

	case 'fatrix2'
		T = fatrix2('arg', arg2, ...
			'imask', arg2.Gmri.mask, ...
			'omask', arg2.Gmri.mask, ...
			'idim', arg2.Gmri.Nd, ...
			'odim', arg2.Gmri.Nd, ...
			'forw', @Gmri_zmap_forw, ...
			'back', @Gmri_zmap_forw); % trick: because Hermitian

	otherwise
		fail('unknown class "%s"', class(arg2.Gmri))
	end
end


% Gmri_gram_new_zmap()
% update Toeplitz-like gram-matrix Fatrix object for new zmap
function T = Gmri_gram_new_zmap(T, varargin) % (ti, zmap, L, aL)
T.arg.Gmri = T.arg.Gmri.arg.new_zmap(T.arg.Gmri, varargin{:}); % yikes!
T = Gmri_gram_work(T.arg, []);


% Gmri_zmap_forw_Fatrix()
% y = T * x
function y = Gmri_zmap_forw_Fatrix(arg2, x)

arg1 = arg2.Gmri.arg;

if size(x,1) ~= arg2.dim(2)
	x = reshapee(x, prod(arg1.Nd), []); % [(N) (nc)] to [*N *nc]
	x = x(arg1.mask,:); % [np *nc]
end
nc = ncol(x);

LL = ncol(arg1.aB);
y = 0;
for ll=1:LL
	tmp = repmat(arg1.aCt(:,ll), [1 nc]) .* x;
	tmp = arg2.T{ll} * tmp;
	tmp = repmat(conj(arg1.aCt(:,ll)), [1 nc]) .* tmp;
	y = y + tmp;
%	y = y + conj(arg1.aCt(:,ll)) .* (arg2.T{ll} * (arg1.aCt(:,ll) .* x));
end


% Gmri_zmap_forw()
function y = Gmri_zmap_forw(arg2, x)
y = Gmri_zmap_forw_Fatrix(arg2, x);
y = embed(y, arg2.Gmri.arg.mask); % required for fatrix2


% Gmri_gram_test()
function Gmri_gram_test
ig = image_geom('nx', 6, 'ny', 8, 'dx', 1, 'offsets', 'dsp');
ig.mask = ellipse_im(ig) > 0;
x = ellipse_im(ig, 'shepplogan-emis', 'type', 'slow');
kspace = mri_trajectory('spiral1', {}, [ig.nx ig.ny], ig.fov);
ti = linspace(0, 10e-3, size(kspace,1));
zmap = 20 * ig.ones + 2i * pi * 10;
L = 8;
%A = Gmri(kspace, ig.mask, 'exact', 1, 'n_shift', ig.dim/2); % perfect gram
%A = Gmri(kspace, ig.mask, 'class', 'fatrix2'); % gram mult not so well matched
A = Gmri(kspace, ig.mask, 'L', L, 'ti', ti, 'zmap', zmap);
T = build_gram(A);
fatrix2_tests(T, 'complex', 1)
x1 = A' * (A * x);
x1 = ig.embed(x1);
x2 = T * x;
im plc 1 3
im(1, x1)
im(2, x2)
im(3, x2 - x1)
try
	equivs(x2, x1, 'thresh', 5e-3) % todo: why so big - due to nufft?
%	equivs(x2, x1)
catch
	keyboard
end
% fatrix2_tests(A, 'complex', 1)
