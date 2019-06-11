 function ob = Gmri_SENSE(kspace, mask, varargin)
%function ob = Gmri(kspace, mask, 'option name', option, ...)
%
% Construct Gmri object for MR image reconstruction, or for other problems
% that require reconstruction of a function from nonuniform Fourier samples.
% This object can also model relaxation or off-resonance effects for MRI.
%
% Essentially, you create a system matrix object by calling:
%	G = Gmri( ... )
% and then you can use it thereafter by typing commands like
%	y = G * x;
% which will auto-magically emulate the MR signal equation.
% This is useful for iterative image reconstruction in MRI.
%
% Besides simple utilities like display, there are the following
% capabilities of this object:
%	y = G * x		forward operation
%	x = G' * y		adjoint operation
%
% in:
%	kspace		[M,d]		trajectory in 1/mm or 1/cm
%					(use same units as 'fov')
%	mask		[[Nd]]		logical support array (see image_geom)
%
% options: (each must be named in quoted string, followed by argument.)
%	'fov'		[d,1]		field of view (typically in mm or cm)
%	'basis'		cell		parameters for image-domain basis:
%						{'rect', [dx dy]} for pixels
%					(default: 'dirac')
%	'nufft'		cell		arguments for nufft_init() via Gnufft()
%	'exact'		1|0		1: use exact (slow) transform.
%					0 (default): use time-seg approximation
%	'n_shift'	[1,d]		image-domain shift (for 'exact' only)
%
% options for field-corrected reconstruction:
%	ti		[d,1]		sample times
%	zmap		[[Nd]]		relax_map + 2i*pi*field_map
%	L		[1]		# of approximation terms
%				or use {Linit, rmsmax} (see mri_exp_approx.m)
%	aL		[1]		# of terms for autocorrelation histogram
%					that is used for Toeplitz version.
%
% After building this object, with or without a zmap, the user can update
% or change the zmap (e.g. for dynamic cases) using the following call:
%	ob = feval(ob.arg.new_zmap, ob, ti, zmap, L, aL);
%
% To change basis "B" and coefficients (transposed) "Ct" for modified models:
%	ob = feval(ob.arg.new_B_Ct, ob, B, Ct);
%	ob = feval(ob.arg.new_aB_aCt, ob, aB, aCt); % for Toeplitz/Gram version
%
% To change image basis (e.g., from dirac for CP to rect for iterative), use:
%	ob = feval(ob.arg.new_image_basis, ob, basis);
% where basis is a cell array of arguments as in 'basis' above.
%
% Copyright 2004-4-22, Jeff Fessler, The University of Michigan

if nargin == 1 & streq(kspace, 'test'), Gmri_test, return, end
if nargin < 5, help(mfilename), error(mfilename), end
if ~islogical(mask), error 'mask must be logical', end

% defaults
arg.segtype = 'hist,time,unif';
arg.exact = 0;
arg.Gnufft = [];
arg.fov = 1;
arg.basis_args = {'dirac'};
arg.basis = [];
arg.zmap = [];
arg.ti = [];
arg.L = [];
arg.aL = []; % for autocorrelation zmap. default is arg.L
arg.sens = [];
N = size(mask);
arg.nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
arg.n_shift = [];
arg.exp_approx_args = {}; % arguments for mri_exp_approx()

for ii=1:length(varargin) % trick: handle old style 'exact' w/o 1/0 argument
	if streq(varargin{ii}, 'exact') & ...
		(ii==length(varargin) | ischar(varargin{ii+1}))
		varargin = {varargin{1:ii}, 1, varargin{(ii+1):end}};
		warning 'old style ''exact'' option: use ''exact'', 1 instead'
		break
	end
end

arg = vararg_pair(arg, varargin, 'subs', ...
	{'basis', 'basis_args'; 'nufft', 'nufft_args'});

if ~isa(arg.basis_args, 'cell'), error 'basis_args must be cell array', end
if ~isa(arg.nufft_args, 'cell'), error 'nufft_args must be cell array', end

if ~isempty(arg.zmap) && ~arg.exact & (isempty(arg.L) | isempty(arg.ti))
	error 'for field-corrected recon, need zmap, ti, and L'
end

if sum(abs(arg.zmap(:))) == 0
  arg.L = [];
  arg.ti = [];
  arg.zmap = [];
end

dd = size(kspace,2);
arg.kspace = kspace;
if length(arg.fov) == 1
	arg.fov = arg.fov * ones(1,dd);
end
arg.mask = mask;	% true(nx,ny);
arg.Nd = size(mask);
arg.dim = [size(kspace,1) sum(mask(:))];

arg.Nc = size(arg.sens,2);

dimp = arg.dim;
if ~isempty(arg.sens)
  dimp = [arg.dim(1)*arg.Nc arg.dim(2)];
end

omega = zeros(size(kspace));
for id=1:dd
	omega(:,id) = 2*pi*kspace(:,id) * arg.fov(id) / arg.Nd(id);
end
if max(abs(omega(:))) > pi+1e-6
	warn 'omega exceeds pi.  are you sure you intended this?'
end

%
% initialize exact transform or NUFFT approximation object
%
if arg.exact
	if isempty(arg.n_shift), error 'n_shift required for exact', end
	if length(N) == 2
		[n1 n2] = ndgrid([0:N(1)-1] - arg.n_shift(1), ...
				[0:N(2)-1] - arg.n_shift(2));
		arg.u = [n1(mask(:)) n2(mask(:))]'; % [2,np]
	elseif length(N) == 3
		[n1 n2 n3] = ndgrid([0:N(1)-1] - arg.n_shift(1), ...
				[0:N(2)-1] - arg.n_shift(2), ...
				[0:N(3)-1] - arg.n_shift(3));
		arg.u = [n1(mask(:)) n2(mask(:)) n3(mask(:))]'; % [3,np]
	else
		error 'not done'
	end
	arg.u = complexify(arg.u);
	arg.v = 1i * omega.';
	[arg.u, arg.v] = single(arg.u, arg.v);
else
	if ~isempty(arg.n_shift), error 'n_shift inapplicable', end
	arg.Gnufft = Gnufft(mask, {omega, arg.nufft_args{:}});
end

% function handle for changing zmap
arg.new_image_basis = @Gmri_new_image_basis;
arg.new_zmap = @Gmri_new_zmap;
arg.new_B_Ct = @Gmri_new_B_Ct;
arg.new_aB_aCt = @Gmri_new_aB_aCt;

%
% build Fatrix object
%
if arg.exact
	ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
		'forw', @Gmri_forw_exact, 'back', @Gmri_back_exact, ...
		'gram', @Gmri_gram);
else
	ob = Fatrix(dimp, arg, 'caller', mfilename, ...
		'forw', @Gmri_forw, 'back', @Gmri_back, 'gram', @Gmri_gram);
end

%
% spatial basis Fourier transform
%
%arg.basis = Gmri_basis(arg.basis_args, arg.dim, arg);
ob = Gmri_new_image_basis(ob, arg.basis_args);

% add zmap effects if available now.
% else the user can use this call to add them later (for dynamic field maps).
% trick: this avoids rebuilding Gnufft for each zmap!
if ~isempty(arg.zmap)
	ob = feval(ob.arg.new_zmap, ob, arg.ti, arg.zmap, arg.L);
end


%
% Gmri_new_image_basis()
%
function ob = Gmri_new_image_basis(ob, basis_args)
arg = ob.arg;
arg.basis = Gmri_basis(basis_args, arg.dim, arg);
if ~isempty(arg.sens)
  arg.basis.transform = repmat(arg.basis.transform,[arg.Nc 1]);
end
ob.arg = arg;


%
% Gmri_basis()
% precompute Fourier transform of spatial basis function
%
function basis = Gmri_basis(basis_args, dim, arg)

basis.args = basis_args;
basis.type = basis_args{1};
dd = size(arg.kspace,2);

switch basis.type
case 'dirac'
	Bi = ones(dim(1), 1);

case 'rect'
	if length(basis_args) == 1
		dx = arg.fov ./ arg.Nd;	% usual default
	else
		if length(basis_args) ~= dd, error 'bad basis_args', end
		dx = basis_args{2};
		if length(dx) == 1, dx = dx * ones(1,dd); end
		if length(dx) ~= dd, error 'bad dx'; end
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

basis.transform = dsingle(Bi);


%
% Gmri_new_zmap()
% update arguments of Gmri object to reflect a new zmap
%
function G = Gmri_new_zmap(G, ti, zmap, L, aL)

if nargin < 3, fail('Gmri_new_zmap needs at least 3 arguments'), end

arg = G.arg;
if length(ti) ~= size(arg.kspace,1), fail('ti size mismatch'), end
if ~isequal(size(zmap), size(arg.mask)), fail('zmap size mismatch'), end
arg.ti = ti(:);
arg.zmap = zmap(arg.mask);

% trick to handle 'exact' case (for which L is irrelevant)
if arg.exact
	dd = length(arg.Nd);
	% trick: if already a zmap in u,v, then replace it
	arg.u = complexify([arg.u(1:dd,:); arg.zmap.']); % [DD+1,np]
	arg.v = complexify([arg.v(1:dd,:); arg.ti']); % [DD+1,nd]

	G.arg = arg;
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

G.arg = arg;


%
% Gmri_new_B_Ct()
% new values of basis "B" and coefficients (transposed) "Ct"
% B [M,L]
% C [L,N] so Ct is [N,L]
%
function G = Gmri_new_B_Ct(G, B, Ct)
G.arg.zmap = true; % trick: make it appear there is a zmap
G.arg.L = ncol(B);
if G.arg.L ~= ncol(Ct), fail('L mismatch'), end
if nrow(B) ~= G.arg.dim(1), fail('data size M mismatch'), end
if nrow(Ct) ~= G.arg.dim(2), fail('image size N mismatch'), end
G.arg.B = B;
G.arg.Ct = Ct;


%
% Gmri_new_aB_aCt()
% new values of basis "aB" and coefficients (transposed) "aCt"
% (for auto-correlation histogram case used in Toeplitz / Gram version)
% aB [M,L]
% aC [L,N] so aCt is [N,L]
%
function G = Gmri_new_aB_aCt(G, aB, aCt)
G.arg.zmap = true; % trick: make it appear there is a zmap
G.arg.L = ncol(aB);
if G.arg.L ~= ncol(Cat), fail('L mismatch'), end
if nrow(aB) ~= G.arg.dim(1), fail('data size M mismatch'), end
if nrow(aCt) ~= G.arg.dim(2), fail('image size N mismatch'), end
G.arg.aB = aB;
G.arg.aCt = aCt;

%
% Gmri_forw_exact(): y = G * x
%
function y = Gmri_forw_exact(arg, x)

if size(x,1) ~= arg.dim(2)
	x = reshape(x, prod(arg.Nd), []); % [(N),(nc)] to [*N,nc]
	x = x(arg.mask,:); % [np,*nc]
end
y = exp_xform_mex(complexify(single(x)), ...
	complexify(single(arg.u)), complexify(single(arg.v)));
y = double(y);
y = y .* repmat(arg.basis.transform, [1 ncol(y)]); % [nd, *nc]


%
% Gmri_forw(): y = G * x
%
function y = Gmri_forw(arg, x)

Nt = arg.dim(1);

if size(x,1) ~= arg.dim(2)
	x = reshape(x, prod(arg.Nd), []); % [(N),(nc)] to [*N,*nc]
	x = x(arg.mask,:); % [np,*nc]
end
nc = ncol(x);

if isempty(arg.zmap)
  if isempty(arg.sens)
    y = arg.Gnufft * x;
  else % parallel imaging
    y = zeros(Nt*arg.Nc,1);
    for cc=1:arg.Nc
      y((cc-1)*Nt+1:cc*Nt) = arg.Gnufft * (arg.sens(:,cc) .* x);
    end
  end
else % approximation
  if isempty(arg.sens)
	y = 0;
	for ll=1:arg.L
		tmp = repmat(arg.Ct(:,ll), [1 nc]) .* x;
		tmp = arg.Gnufft * tmp;
		tmp = repmat(arg.B(:,ll), [1 nc]) .* tmp;
		y = y + tmp;
%		y = y + arg.B(:,ll) .* (arg.Gnufft * (arg.Ct(:,ll) .* x));
	end
  else % parallel imaging
    y = zeros(Nt*arg.Nc,1);
    for ll=1:arg.L
      for cc=1:arg.Nc
        tmp = repmat(arg.sens(:,cc).*arg.Ct(:,ll), [1 nc]) .*x;
        tmp = arg.Gnufft * tmp;
        tmp = repmat(arg.B(:,ll), [1 nc]) .* tmp;
        y((cc-1)*Nt+1:cc*Nt) = y((cc-1)*Nt+1:cc*Nt) + tmp;
      end
    end
  end
end
y = y .* double(repmat(arg.basis.transform, [1 ncol(y)])); % [nd, *nc]


%
% Gmri_back_exact(): x = G' * y
%
function x = Gmri_back_exact(arg, y)
y = y .* repmat(conj(arg.basis.transform), [1 ncol(y)]);

% trick: conj(exp(-uv)) = exp(-conj(u) conj(v))
vc = complexify(single(conj(arg.v)));
uc = complexify(single(conj(arg.u)));
x = exp_xform_mex(complexify(single(y)), vc, uc);
x = double(x);


%
% Gmri_back(): x = G' * y
% full adjoint ("back-projection")
%
function x = Gmri_back(arg, y)
nc = ncol(y);
Nt = arg.dim(1);
y = y .* repmat(conj(arg.basis.transform), [1 nc]);

if isempty(arg.zmap)
  if isempty(arg.sens)
	x = arg.Gnufft' * y;
  else % parallel imaging
    x = 0;
    for cc=1:arg.Nc
      tmp = y((cc-1)*Nt+1:cc*Nt);
      if any(tmp)
        x = x + conj(arg.sens(:,cc)) .* (arg.Gnufft' * tmp);
      end
    end
  end
else % approximation
  if isempty(arg.sens)
	x = 0;
	for ll=1:arg.L
		tmp = repmat(conj(arg.B(:,ll)), [1 nc]) .* y;
		tmp = arg.Gnufft' * tmp;
		tmp = repmat(conj(arg.Ct(:,ll)), [1 nc]) .* tmp;
		x = x + tmp;
%		x = x + conj(arg.Ct(:,ll)) .* (arg.Gnufft' * (conj(arg.B(:,ll)) .* y));
	end
  else % parallel imaging
    x = 0;
    for ll=1:arg.L
      for cc=1:arg.Nc
        tmp = y((cc-1)*Nt+1:cc*Nt);
        if any(tmp)
          tmp = repmat(conj(arg.B(:,ll)), [1 nc]) .* tmp;
          tmp = arg.Gnufft' * tmp;
          tmp = repmat(conj(arg.sens(:,cc)).*conj(arg.Ct(:,ll)), [1 nc]) .* tmp;
          x = x + tmp;
        end
      end
    end
  end
end


%
% Gmri_test
%
function Gmri_test

% test cases with tiny objects
ig = image_geom('nx', 8, 'ny', 6, 'dx', 1, 'offsets', 'dsp');
ig.mask = ellipse_im(ig) > 0;
kspace = mri_trajectory('spiral1', {}, [ig.nx ig.ny], ig.fov);
ti = linspace(0, 10e-3, size(kspace,1));
zmap = 20 * ig.ones + 2i * pi * 10;

% exact
if 1
	G = Gmri(kspace, ig.mask, 'exact', true, ...
        	'ti', ti, 'zmap', zmap, 'n_shift', [0 0]);
	G = feval(G.arg.new_image_basis, G, {'dirac'});
	
	Fatrix_test_basic(G, ig.mask)
end

% usual
if 1
	G = Gmri(kspace, ig.mask, 'L', 8, ...
        	'ti', ti, 'zmap', zmap);

	Fatrix_test_basic(G, ig.mask)
end

% toeplitz
if 1
	G = G.arg.new_zmap(G, ti, ig.ones * 2i * pi * 10, 4);
	T = build_gram(G, 1);
	Fatrix_test_basic(T, ig.mask)

	if 1
		t1 = G' * G; t1 = t1(:,:);
		t2 = T(:,:);
	 	max_percent_diff(t1,t2)
	end
	T = T.arg.new_zmap(T, ti, ig.ones * 2i * pi * 5, 4); % test
end

% other external tests
Gmri_test_exact

printm('to further test this, running mri_example.m')
run_mfile_local('mri_example')
