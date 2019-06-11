 function ob = Gmri(kspace, mask, varargin)
%function ob = Gmri(kspace, mask, 'option name', option, ...)
%|
%| Construct Gmri object for MR image reconstruction, or for other problems
%| that require reconstruction of a function from nonuniform Fourier samples.
%| This object can also model relaxation and/or off-resonance effects for MRI.
%|
%| Essentially, you create a system matrix object by calling:
%|	A = Gmri( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically emulate the MR signal equation.
%| This is useful for iterative image reconstruction in MRI.
%|
%| Besides simple utilities like display, there are the following
%| capabilities of this object:
%|	y = A * x		forward operation
%|	x = A' * y		adjoint operation
%|
%| in
%|	kspace		[M d]		trajectory in 1/mm or 1/cm
%|					(use same units as 'fov')
%|	mask		[[Nd]]		logical support array (see image_geom)
%|
%| options
%|	'fov'		[d 1]		field of view (typically in mm or cm)
%|					default: 1
%|	'basis'		cell		parameters for image-domain basis:
%|						{'rect', [dx dy]} for pixels
%|					(default: {'dirac'})
%|	'nufft'		cell		arguments for nufft_init() via Gnufft()
%|					(excluding the first argment, omega)
%|	'exact'		1|0		1: use exact (slow) transform.
%|					0 (default): use default in mri_exp_approx.m
%|	'n_shift'	[1 d]		image-domain shift (for 'exact' only)
%|						otherwise put in 'nufft' args!
%|
%| options for field-corrected reconstruction:
%|	ti		[d 1]		sample times
%|	zmap		[[Nd]]		relax_map + 2i*pi*field_map
%|					included as a exp(-zmap * ti) term!
%|	L		[1]		# of approximation terms
%|				or use {Linit, rmsmax} (see mri_exp_approx.m)
%|	aL		[1]		# of terms for autocorrelation histogram
%|					that is used for Toeplitz version.
%|	exp_approx_args	{}		arguments for mri_exp_approx()
%|
%| out
%|	ob		[nd np]		fatrix2 or Fatrix object
%|
%| After building this object, with or without a zmap, the user can update
%| or change the zmap (e.g. for dynamic cases) using the following call:
%|	ob = ob.arg.new_zmap(ob, ti, zmap, L, aL);
%|
%| To change basis "B" and coefficients (transposed) "Ct" for modified models:
%|	ob = ob.arg.new_B_Ct(ob, B, Ct);
%|	ob = ob.arg.new_aB_aCt(ob, aB, aCt); % for Toeplitz/Gram version
%|
%| To change image basis (e.g., from dirac for CP to rect for iterative), use:
%|	ob = ob.arg.new_image_basis(ob, basis);
%| where basis is a cell array of arguments as in 'basis' above.
%|
%| Copyright 2004-4-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(kspace, 'test'), Gmri_test, return, end
if nargin < 2, ir_usage, end
if ~islogical(mask), fail 'mask must be logical', end

% defaults
arg.class = 'fatrix2';
%arg.class = 'Fatrix';
arg.exact = 0;
arg.Gnufft = [];
arg.fov = 1;
arg.basis_args = {'dirac'};
arg.basis = [];
arg.zmap = [];
arg.ti = [];
arg.L = [];
arg.aL = []; % for autocorrelation zmap. default is arg.L
N = size(mask);
arg.nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
arg.n_shift = [];
arg.exp_approx_args = {}; % arguments for mri_exp_approx()

for ii=1:length(varargin) % trick: handle old style 'exact' w/o 1/0 argument
	if streq(varargin{ii}, 'exact') && ...
		(ii==length(varargin) || ischar(varargin{ii+1}))
		varargin = {varargin{1:ii}, 1, varargin{(ii+1):end}};
		warning 'old style ''exact'' option: use ''exact'', 1 instead'
		break
	end
end

arg = vararg_pair(arg, varargin, 'subs', ...
	{'basis', 'basis_args'; 'nufft', 'nufft_args'});

if ~isa(arg.basis_args, 'cell'), fail 'basis_args must be cell array', end
if ~isa(arg.nufft_args, 'cell'), fail 'nufft_args must be cell array', end

if ~isempty(arg.zmap) && ~arg.exact && (isempty(arg.L) || isempty(arg.ti))
	fail 'for field-corrected recon, need zmap, ti, and L'
end

dd = size(kspace,2);
arg.kspace = kspace;
if length(arg.fov) == 1
	arg.fov = arg.fov * ones(1,dd);
end
arg.mask = mask;
arg.Nd = size(mask);
arg.dim = [size(kspace,1) sum(mask(:))];

omega = zeros(size(kspace), class(kspace));
for id=1:dd
	omega(:,id) = 2*pi*kspace(:,id) * arg.fov(id) / arg.Nd(id);
end
if max(abs(omega(:))) > pi+1e-6
	warn 'omega exceeds pi.  are you sure you intended this?'
end

% initialize exact transform or NUFFT approximation object
if arg.exact
	if isempty(arg.n_shift), fail 'n_shift required for exact', end
	if length(N) == 2
		if N(2) == 1
			if length(arg.n_shift) ~= 1, fail 'n_shift size', end
			n1 = [0:N(1)-1]' - arg.n_shift(1);
			arg.u = n1(mask(:))'; % [1 np]
		else
			if length(arg.n_shift) ~= 2, fail 'n_shift size', end
			[n1 n2] = ndgrid([0:N(1)-1] - arg.n_shift(1), ...
					[0:N(2)-1] - arg.n_shift(2));
			arg.u = [n1(mask(:)) n2(mask(:))]'; % [2 np]
		end
	elseif length(N) == 3
		if length(arg.n_shift) ~= 3, fail 'n_shift size', end
		[n1 n2 n3] = ndgrid([0:N(1)-1] - arg.n_shift(1), ...
				[0:N(2)-1] - arg.n_shift(2), ...
				[0:N(3)-1] - arg.n_shift(3));
		arg.u = [n1(mask(:)) n2(mask(:)) n3(mask(:))]'; % [3 np]
	else
		fail 'not done'
	end
	arg.v = 1i * omega.';
	arg.u = single(arg.u);
	arg.u = complexify(arg.u);
	arg.v = single(arg.v);
else
	if ~isempty(arg.n_shift), fail 'n_shift inapplicable', end
	arg.Gnufft = Gnufft(mask, {omega, arg.nufft_args{:}});
end

% function handle for changing zmap
arg.new_image_basis = @Gmri_new_image_basis;
arg.new_zmap = @Gmri_new_zmap;
arg.new_B_Ct = @Gmri_new_B_Ct;
arg.new_aB_aCt = @Gmri_new_aB_aCt;

% build object
switch arg.class
case 'Fatrix'
	if arg.exact
		ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
			'forw', @Gmri_forw_exact_Fatrix, ...
			'back', @Gmri_back_exact_Fatrix, ...
			'gram', @Gmri_gram);
	else
		ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
			'forw', @Gmri_forw, ...
			'back', @Gmri_back_Fatrix, ...
			'gram', @Gmri_gram);
	end

case 'fatrix2'
	if arg.exact
		ob = fatrix2('mask', mask, 'arg', arg, ...
			'odim', arg.dim(1), 'does_many', true, ...
			'forw', @Gmri_forw_exact, 'back', @Gmri_back_exact);
	else
		ob = fatrix2('mask', mask, 'arg', arg, ...
			'odim', arg.dim(1), 'does_many', true, ...
			'forw', @Gmri_forw, 'back', @Gmri_back, ...
			'gram', @Gmri_gram);
	end

otherwise
	fail('class "%s" unknown', arg.class)
end


% spatial basis Fourier transform
%arg.basis = Gmri_basis(arg.basis_args, arg.dim, arg);
ob = Gmri_new_image_basis(ob, arg.basis_args);

% add zmap effects if available now.
% else the user can use this call to add them later (for dynamic field maps).
% trick: this avoids rebuilding Gnufft for each zmap!
if ~isempty(arg.zmap)
	ob = ob.arg.new_zmap(ob, arg.ti, arg.zmap, arg.L);
end


% Gmri_new_image_basis()
function ob = Gmri_new_image_basis(ob, basis_args)
arg = ob.arg;
ob.arg.basis = Gmri_basis(basis_args, arg.dim, arg);
%arg.basis = Gmri_basis(basis_args, arg.dim, arg);
%ob.arg = arg;


% Gmri_basis()
% precompute Fourier transform of spatial basis function
function basis = Gmri_basis(basis_args, dim, arg)

basis.args = basis_args;
basis.type = basis_args{1};
dd = size(arg.kspace,2);

switch basis.type
case 'dirac'
	Bi = ones(dim(1), 1, 'single');

% 'sinc' and 'dirac*dx' simply provide an appropriate "scale factor" dx
% to relate Fourier integral and summation
case {'sinc', 'dirac*dx'}
	Bi = ones(dim(1), 1, 'single') * prod(arg.fov ./ arg.Nd);

case 'rect' % rect(x/dx) <=> |dx| * sinc(dx * u)
	if length(basis_args) == 1
		dx = abs(arg.fov) ./ arg.Nd; % usual default
	else
		if length(basis_args) ~= dd, fail 'bad basis_args', end
		dx = abs(basis_args{2});
		if length(dx) == 1, dx = dx * ones(1,dd); end
		if length(dx) ~= dd, fail 'bad dx'; end
	end
	Bi = ones(dim(1),1);
	for id=1:size(arg.kspace,2)
		if dx(id) % allow for zero-sized pixels (Dirac impulses)
			Bi = Bi .* (dx(id) * nufft_sinc(dx(id) * arg.kspace(:,id)));
		end
	end

otherwise
	fail('unknown basis_type "%s"', arg.basis_type)
end

basis.transform = single(Bi);


% Gmri_new_zmap()
% update arguments of Gmri object to reflect a new zmap
% A = Gmri_new_zmap(A, ti, zmap, L, aL)
function A = Gmri_new_zmap(A, ti, zmap, varargin)
if nargin < 3, fail('Gmri_new_zmap needs at least 3 arguments'), end
arg = Gmri_new_zmap_arg(A.arg, ti, zmap, varargin{:});
fields = {'ti', 'zmap', 'u', 'v', 'B', 'Ct', 'L', 'aB', 'aCt', 'aL'};
for ii=1:numel(fields)
	tmp = fields{ii};
	if isfield(arg, tmp)
		A.arg.(tmp) = arg.(tmp);
	end
end


% Gmri_new_zmap_arg()
% update arguments of Gmri object to reflect a new zmap
% side effects arg.{ti zmap u v B Ct L aB aCt aL}
function arg = Gmri_new_zmap_arg(arg, ti, zmap, L, aL)

if length(ti) ~= size(arg.kspace,1), fail('ti size mismatch'), end
if ~isequal(size(zmap), size(arg.mask)), fail('zmap size mismatch'), end
arg.ti = ti(:);
arg.zmap = zmap(arg.mask);

% trick to handle 'exact' case (for which L is irrelevant)
if arg.exact
	dd = length(arg.Nd);
	if arg.Nd(2) == 1 % 1d
		dd = 1;
	end
	% trick: if already a zmap in u,v, then replace it
	arg.u = complexify([arg.u(1:dd,:); arg.zmap.']); % [DD+1 np]
	arg.v = complexify([arg.v(1:dd,:); arg.ti']); % [DD+1 nd]
return
end

if ~isvar('L') || isempty(L) % 4th argument is optional, so defaults here:
	if ~isempty(arg.L)
		L = arg.L;
	else
		fail('user must provide arg.L or L input')
	end
end

% initialize exponential approximations for field-corrected reconstruction
if any(real(arg.zmap))
	arg.exp_approx_type = {'hist,time,unif', [40 10]};
else
	arg.exp_approx_type = {'hist,time,unif', 40};
end
[arg.B C] = mri_exp_approx(arg.ti, arg.zmap, L, ...
		'type', arg.exp_approx_type, arg.exp_approx_args{:});
if any(isnan(arg.B(:)))
	warning 'bug: nan values in B'
	keyboard
end
arg.Ct = C.';

if ~isvar('aL') || isempty(aL) % 5th argument is optional, so defaults here:
	if ~isempty(arg.aL)
		aL = arg.aL;
	elseif iscell(L) % trick: use "found" L because aL >= L generally
		aL = {ncol(arg.B), L{2}}; % {L, rmsmax}
	else
		aL = ncol(arg.B);
	end
end

% store size (seems redundant, but OK)
arg.L = ncol(arg.B);
if iscell(L)
	arg.rmsmax = L{2};
	printm('L=%d found', arg.L)
end

% generate one with auto-correlation histogram too.
% this may not be needed, but the time is small compared to Gnufft
if any(real(arg.zmap) ~= 0) && all(imag(arg.zmap) == 0)
	warn('zmap has only a real component!? are you sure?')
end

[arg.aB C] = mri_exp_approx(arg.ti, arg.zmap, aL, 'acorr', 1, ...
	'type', arg.exp_approx_type, arg.exp_approx_args{:});

% trick: for symmetric histogram, it should be real!
arg.aB = reale(arg.aB, 'warn', 'in Gmri');
arg.aCt = C.';

arg.aL = ncol(arg.aB);
if iscell(aL)
	printm('aL=%d found', arg.aL)
end


% Gmri_new_B_Ct()
% new values of basis "B" and coefficients (transposed) "Ct"
% B [M L]
% C [L N] so Ct is [N L]
%
function A = Gmri_new_B_Ct(A, B, Ct)
A.arg.zmap = true; % trick: make it appear there is a zmap
A.arg.L = ncol(B);
if A.arg.L ~= ncol(Ct), fail('L mismatch'), end
if nrow(B) ~= A.arg.dim(1), fail('data size M mismatch'), end
if nrow(Ct) ~= A.arg.dim(2), fail('image size N mismatch'), end
A.arg.B = B;
A.arg.Ct = Ct;


% Gmri_new_aB_aCt()
% new values of basis "aB" and coefficients (transposed) "aCt"
% (for auto-correlation histogram case used in Toeplitz / Gram version)
% aB [M L]
% aC [L N] so aCt is [N L]
%
function A = Gmri_new_aB_aCt(A, aB, aCt)
A.arg.zmap = true; % trick: make it appear there is a zmap
A.arg.L = ncol(aB);
if A.arg.L ~= ncol(Cat), fail('L mismatch'), end
if nrow(aB) ~= A.arg.dim(1), fail('data size M mismatch'), end
if nrow(aCt) ~= A.arg.dim(2), fail('image size N mismatch'), end
A.arg.aB = aB;
A.arg.aCt = aCt;


% Gmri_forw_exact_Fatrix(): y = A * x
function y = Gmri_forw_exact_Fatrix(arg, x)

if size(x,1) ~= arg.dim(2)
	x = reshape(x, prod(arg.Nd), []); % [(N) (nc)] to [*N nc]
	x = x(arg.mask,:); % [np *nc]
end
y = exp_xform_mex(complexify(single(x)), ...
	complexify(arg.u), complexify(arg.v));
y = y .* repmat(arg.basis.transform, [1 ncol(y)]); % [nd *nc]


% Gmri_forw_exact(): y = A * x
function y = Gmri_forw_exact(arg, x)
x = reshapee(x, prod(arg.Nd), []); % [(N) (nc)] to [*N *nc]
x = x(arg.mask,:); % [np *nc]
y = exp_xform_mex(complexify(single(x)), arg.u, arg.v); % [nd *nc]
y = y .* repmat(arg.basis.transform, [1 ncol(y)]); % [nd *nc]


% Gmri_forw(): y = A * x
function y = Gmri_forw(arg, x)

if size(x,1) ~= arg.dim(2)
	x = reshape(x, prod(arg.Nd), []); % [(N) (nc)] to [*N *nc]
	x = x(arg.mask,:); % [np *nc]
end
nc = ncol(x);

if isempty(arg.zmap)
	y = arg.Gnufft * x;
else % approximation
	y = 0;
	for ll=1:arg.L
		tmp = repmat(arg.Ct(:,ll), [1 nc]) .* x;
		tmp = arg.Gnufft * tmp;
		tmp = repmat(arg.B(:,ll), [1 nc]) .* tmp;
		y = y + tmp;
%		y = y + arg.B(:,ll) .* (arg.Gnufft * (arg.Ct(:,ll) .* x));
	end
end
y = y .* repmat(arg.basis.transform, [1 nc]); % [nd *nc]


% Gmri_back_exact_Fatrix(): x = A' * y
function x = Gmri_back_exact_Fatrix(arg, y)
y = y .* repmat(conj(arg.basis.transform), [1 ncol(y)]); % [nd nc]
% trick: conj(exp(-uv)) = exp(-conj(u) conj(v))
vc = complexify(conj(arg.v));
uc = complexify(conj(arg.u));
x = exp_xform_mex(complexify(single(y)), vc, uc); % [np nc]


% Gmri_back_exact(): x = A' * y
function x = Gmri_back_exact(arg, y)
x = Gmri_back_exact_Fatrix(arg, y);
x = embed(x, arg.mask); % [(N) nc] as required for fatrix2


% Gmri_back_Fatrix(): x = A' * y
% full adjoint ("back-projection")
function x = Gmri_back_Fatrix(arg, y)
nc = ncol(y);
y = y .* repmat(conj(arg.basis.transform), [1 nc]);

if isempty(arg.zmap)
	x = arg.Gnufft' * y;
else % approximation
	x = 0;
	for ll=1:arg.L
		tmp = repmat(conj(arg.B(:,ll)), [1 nc]) .* y;
		tmp = arg.Gnufft' * tmp;
		tmp = repmat(conj(arg.Ct(:,ll)), [1 nc]) .* tmp;
		x = x + tmp;
	end
end


% Gmri_back(): x = A' * y
% full adjoint ("back-projection")
function x = Gmri_back(arg, y)
x = Gmri_back_Fatrix(arg, y);
x = embed(x, arg.mask); % [(N) nc] as required for fatrix2


% Gmri_test()
function Gmri_test

% test cases with tiny objects
ig = image_geom('nx', 8, 'ny', 6, 'dx', 1, 'offsets', 'dsp');
ig.mask = ellipse_im(ig) > 0;
kspace = mri_trajectory('spiral1', {}, [ig.nx ig.ny], ig.fov);
ti = linspace(0, 10e-3, size(kspace,1));
zmap = 20 * ig.ones + 2i * pi * 10;

classes = {'fatrix2', 'Fatrix'};
for ic=1:2
	args = {kspace, ig.mask, 'class', classes{ic}};

	% exact
	if 1
		A = Gmri(args{:}, 'exact', true, ...
			'ti', ti, 'zmap', zmap, 'n_shift', [0 0]);
		A = A.arg.new_image_basis(A, {'dirac'});
		fatrix2_tests(A, 'complex', 1)
	end

	% usual
	if 1
		A = Gmri(args{:}, 'L', 8, 'ti', ti, 'zmap', zmap);
		fatrix2_tests(A, 'complex', 1, 'tol_gram', 5e-3) % trick
	end

	% toeplitz
	if 1
		A = A.arg.new_zmap(A, ti, ig.ones * 2i * pi * 10, 4);
		T = build_gram(A, 1);
		fatrix2_tests(T, 'complex', 1, 'full', 0)

		if 1
			t1 = A' * A; t1 = t1(:,:);
			t2 = T(:,:);
			max_percent_diff(t1,t2)
		end
		T = T.arg.new_zmap(T, ti, ig.ones * 2i * pi * 5, 4); % test
	end
end

% other external tests
Gmri_test_exact

printm('to further test this, running mri_example.m')
run_mfile_local('mri_example')
