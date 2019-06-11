 function st = newfft(om_in, Nd_in, varargin)
%function st = newfft(om, Nd, [options])
%|
%| New version of NUFFT (pun intended) that uses real interpolation kernels.
%| (The original NUFFT code used complex interpolation needlessly.)
%|
%| This returns a "strum" object with methods for both forward and adjoint
%| d-dimensional NUFFT operations.
%| The forward operation is:
%| X(om_m) = \sum_{n=0}^{N-1} x[n] exp(-1i * om_m * n) for m=1,...,M
%| The adjoint operation is:
%| x_adj[n] = \sum_{m=1}^M X(om_m) exp(+1i * om_m * n) for n=0,...,N-1
%| Note that the adjoint is not the "inverse" NUFFT in general.
%|
%| This routine has numerous options for investigative purposes, but is
%| designed so that the default options should be very good choices.
%| Providing the frequencies and the image size and should suffice.
%| Reducing the neighborhood size 'Jd' (default 6) and/or reducing the
%| over-sampled DFT size 'Kd' (default 2*Nd) may be useful for acceleration.
%|
%| in
%|	om [M,d]	"digital" frequencies in radians (can be empty!)
%|			(if empty, then user must provide 'om' to methods)
%|	Nd [d]		image dimensions (N1,N2,...,Nd)
%|
%| options
%|	'Jd' [d]	# of neighbors used (in each direction). (default: 6)
%|	'Kd' [d]	FFT sizes (should be >= Nd). (default: 2*Nd)
%|	n_shift [d]	n = 0-n_shift to N-1-n_shift (default: 0)
%|		Like fft(), the NUFFT expects the signals to be x(0,0), ...
%|		Use n_shift = [N1/2, N2/2, ...] for x(-N1/2,-N2/2,...), ...
%|
%|	'mode'	char	how to compute the NUFFT (default: 'table1')
%|			'exact'	slow FT
%|			'table0' table with nearest neighbor interpolation
%|			'table1' table with linear interpolation (default)
%|			'gg' gaussian factorization of Greengard & Lee 04 (todo)
%|			'sparse' precompute large sparse matrix
%|				caution: this option may require lots of memory!
%|				(this option used internally to make tables too)
%|			recommended: table0 or table1 to save memory.
%|
%|	'oversample' int table oversampling factor (default: 2^9 for table1)
%|	'gram' 0|1	if 1, precompute additional terms needed for gram matrix
%|	'printmem' 0|1	if 1, print memory usage
%|	'phasing' char	'real' : new real kernels (strongly recommended default)
%|			'complex' : original nufft.m complex kernels
%|			'none' : no phase (straw man; table uses it internally)
%|			'flipreal' : real kernels with sign flips (unsupported)
%|	'dotime' char	report cpu time? (default: '')
%|
%|	'ktype'	char	type of interpolation kernel (default: 'minmax:kb')
%|
%| ktype options:
%|	'diric'		Dirichlet interpolator (exact only if Jd = Kd)
%|	'linear'	linear interpolator (a terrible straw man)
%|	'minmax:kb'	minmax interpolator with excellent KB scaling!
%|	'minmax:tuned'	minmax interpolator, somewhat numerically tuned scaling
%|	'minmax:unif'	minmax with uniform scaling factors (not recommended)
%|	'minmax:user'	minmax interpolator with user parameters required:
%|				'alpha', {alpha}, 'beta', {beta}
%|	'kb:minmax'	kaiser-bessel (KB) interpolator (minmax best alpha, m)
%|	'kb:beatty'	KB with parameters from Beatty et al T-MI Jun 2005
%|	'kb:user'	KB with user-specified KB parameters required:
%|				'kb_m', [m] 'kb_alf', [alpha]
%|	@kernel		user-provided inline interpolation kernel(k,J)
%|			(or a cell array of kernels, one for each dimension)
%|			example ..., 'table', 2^11, 'minmax:kb'
%|			where 2^11 is the table over-sampling factor.
%|
%| out
%|	st	strum object with several methods including the following:
%|		st.fft(x, [om])
%|		st.adj(Xo, [om])
%|		st.p		[M, *Kd]	sparse interpolation matrix
%|						(or empty if table-based)
%|		st.sn		[(Nd)]		scaling factors
%|		st.Nd,Jd,Kd,om	copies of inputs
%|
%| *Nd is shorthand for prod(Nd).
%| (Nd) is shorthand for (N1,N2,...,Nd)
%|
%| Copyright 2007-6-3, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(om_in, 'test'), newfft_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

cpu etic

% inputs
st.om = om_in; % [M,d] frequency samples
st.Nd = Nd_in; % [d] signal dimentions
st.dd = length(st.Nd); % dimensionality of input space (usually 2 or 3)

% defaults
st.gram = false;
st.Jd = 6 * ones(1, st.dd);
st.Kd = 2 * st.Nd;
st.n_shift = zeros(1, st.dd);
st.mode = 'table1';
st.oversample = []; % aka Ld for table mode
st.order = []; % for table mode
st.ktype = 'kb:minmax';
st.kb_m = []; % [dd] KB parameters
st.kb_alf = [];
st.alpha = {}; % [dd] minmax parameters
st.beta = {};
st.tol = 0;
st.printmem = false;
st.is_kaiser_scale = false;
st.phasing = 'real'; % use new real table by default
st.dotime = '';

% options
st = vararg_pair(st, varargin);
dotime = st.dotime; st = rmfield(st, 'dotime');

st.phase_before = []; % place holders
st.phase_after = [];
st.flips = [];

% special cases of input sampling pattern
if ischar(st.om)
	st.om = nufft_samples(st.om, st.Nd);
end

% checks
if st.dd ~= length(st.Jd) || st.dd ~= length(st.Kd)
	error 'inconsistent dim'
end
if st.dd ~= length(st.n_shift)
	fail('n_shift needs %d columns', st.dd)
end

if ~isempty(st.om) && st.dd ~= size(st.om,2)
	fail('omega needs %d columns', st.dd)
end

if st.gram, fail 'todo: gram not done', end

if any(st.Kd < st.Nd), warning 'Kd < Nd unlikely to work.  Try Kd=2*Nd', end

%
% "midpoint" of scaling factors
%
switch st.phasing
case 'real'
	st.Nmid = floor(st.Nd / 2); % new
otherwise
	st.Nmid = (st.Nd - 1) / 2; % old
end

%
% different interpolation modes
%
switch st.mode

case 'exact' % exact interpolation for testing
	st = strum(st, { ...
		'fft', @newfft_exact_for, '(x, [om])';
		'adj', @newfft_exact_adj, '(X, [om])';
		'p', @newfft_exact_p, '([om])';
		'sn', @(st) ones(st.Nd), '()';
		});

case {'table0', 'table1'} % precomuted interpolator table
	st = newfft_table_init(st);

	% set up phase corrections not included in table initialization
	if streq(st.phasing, 'real') || streq(st.phasing, 'flipreal')
		st.phase_before = newfft_phase_before(st.Kd, st.Nmid);
		st.phase_after = @(om) newfft_phase_after(om, st.Nmid, st.n_shift);
		if ~isempty(st.om) % phase that goes after interpolation:
			st.phase_after = st.phase_after(st.om);
		end
	end

case 'gg' % greengard's gaussian
	error 'todo: gg'
	st = strum(st, { ...
		'fft', @newfft_gg_for, '(x, [om])';
		'adj', @newfft_gg_adj, '(X, [om])';
		'p', @newfft_gg_p, '([om])';
		'sn', @newfft_gg_sn, '()';
		});

case 'sparse' % interpolator based on a sparse matrix
	if isempty(st.om), error 'sparse mode requires "om"', end
	st = newfft_init_sparse(st);

otherwise
	fail('unknown mode %s', st.mode)
end

if ~isempty(dotime)
	tmp = whos('st');
	tmp = num2str(tmp.bytes);
	tmp = ['newfft setup ' st.mode ' ' st.phasing(1) ' ' tmp ' ' dotime];
	cpu('etoc', tmp)
end


%
% newfft_init_sparse()
%
% create an interpolator based on a sparse matrix.
% this matrix will be large for intersting problem sizes, so this
% mode is not recommended.  but it is supported in part because it
% is needed for generating samples of the interpolator for the table mode.
%
function st = newfft_init_sparse(st)

om = st.om;

%
% different interpolation kernel mechanisms
%
ktype = st.ktype;

switch class(ktype)
case 'cell' % cell array of kernel functions: {kernel1, kernel2, ..., kernelD}
	if isa(ktype{1}, 'inline') || isa(ktype{1}, 'function_handle')
		if length(ktype) ~= dd, error 'wrong # of kernels', end
		st.kernel = ktype;
		ktype = 'inline';
	else
		error 'cell array should be inline kernels!?'
	end

case {'inline', 'function_handle'} % single inline kernel for all dimensions
	for id = 1:st.dd
		st.kernel{id} = ktype; % all same
	end
	ktype = 'inline';

case 'char'
	% a string that describes the type of interpolator, see below

otherwise
	fail('unknown kernel type %s', class(ktype))
end

%
% interpolator set up
%
Nd = st.Nd;
Jd = st.Jd;
Kd = st.Kd;
dd = st.dd;
is_kaiser_scale = false;

switch ktype
case 'inline'
	% already did it above

case 'diric' % exact interpolator
	if any(Jd ~= Kd), warn 'diric inexact unless Jd=Kd', end
	ktype = 'inline';
	for id = 1:dd
		N = Nd(id);
		K = Kd(id);
		if 1 && streq(st.phasing, 'real')
			N = 2 * floor((K+1)/2) - 1; % trick
		end
		st.kernel{id} = @(k,J) N / K * nufft_diric(k, N, K, true);
	end

case 'linear' % linear interpolator straw man
	ktype = 'inline';
	kernel = inline('(1 - abs(k/(J/2))) .* (abs(k) < J/2)', 'k', 'J');
	for id = 1:dd
		st.kernel{id} = kernel;
	end

case 'kb:beatty' % KB with Beatty et al parameters
	is_kaiser_scale = true;
	if ~isempty(st.kb_alf) || ~isempty(st.kb_m)
		warn 'kb_alf and kb_m ignored'
	end
	K_N = Kd ./ Nd;
	st.kb_alf = pi * sqrt( Jd.^2 ./ K_N.^2 .* (K_N - 1/2).^2 - 0.8 );
%	pr st.kb_alf ./ Jd % approximately 2.34 for K_N = 2
	st.kb_m = zeros(1,dd);
	for id = 1:dd
		st.kernel{id} = kaiser_bessel('inline', Jd(id), ...
				st.kb_alf(id), st.kb_m(id));
	end

case 'kb:minmax' % KB with minmax-optimized parameters
	is_kaiser_scale = true;

	if ~isempty(st.kb_alf) || ~isempty(st.kb_m)
		warn 'kb_alf and kb_m ignored'
	end
	for id = 1:dd
		[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
			kaiser_bessel('inline', Jd(id));
	end

case 'kb:user' % KB with user-defined parameters
	is_kaiser_scale = true;

	if isempty(st.kb_alf) || isempty(st.kb_m)
		fail 'kb_alf and kb_m required'
	end
	if (length(st.kb_alf) ~= dd) || (length(st.kb_m) ~= dd)
		fail('#alpha=%d #m=%d vs dd=%d', ...
			length(st.kb_alf), length(st.kb_m), dd)
	end
	for id = 1:dd
		st.kernel{id} = kaiser_bessel('inline', Jd(id), ...
				st.kb_alf(id), st.kb_m(id));
	end

case 'minmax:kb' % minmax interpolator with KB scaling factors
	for id = 1:dd
		[st.alpha{id}, st.beta{id}] = ...
			nufft_alpha_kb_fit(Nd(id), Jd(id), Kd(id), ...
				'Nmid', st.Nmid(id));
	end

case 'minmax:tuned' % minmax with numerically "tuned" scaling factors
	for id = 1:dd
		[st.alpha{id}, st.beta{id}, ok] = ...
			nufft_best_alpha(Jd(id), 0, Kd(id)/Nd(id));
		if ~ok, error 'unknown J,K/N', end
	end

case 'minmax:user' % minmax interpolator with user-provided scaling factors
	if isempty(st.alpha) || isempty(st.beta)
		error 'user must provide alpha/beta'
	end
	if length(st.alpha) ~= dd || length(st.beta) ~= dd
		error 'alpha/beta size mismatch'
	end

case 'minmax:unif' % minmax with straw man uniform scaling factors
	for id = 1:dd
		st.alpha{id} = 1;
		st.beta{id} = 0;
	end

otherwise
	fail('unknown kernel type %s', ktype)
end

%
% scaling factors: "outer product" of 1D vectors
%
st.sn = 1;
for id=1:dd
	if 1 && streq(st.ktype, 'linear')
		tmp = newfft_scale_tri(Nd(id), Jd(id), Kd(id), st.Nmid);
	elseif streq(st.ktype, 'diric')
		tmp = ones(Nd(id),1);
	elseif is_kaiser_scale
		nc = [0:Nd(id)-1]' - st.Nmid(id);
		tmp = 1 ./ kaiser_bessel_ft(...
			nc/Kd(id), Jd(id), st.kb_alf(id), st.kb_m(id), 1);
	elseif streq(ktype, 'inline')
		tmp = 1 ./ nufft_interp_zn(0, Nd(id), Jd(id), Kd(id), ...
			st.kernel{id}, st.Nmid(id));
	else
		tmp = nufft_scale(Nd(id), Kd(id), ...
			st.alpha{id}, st.beta{id}, st.Nmid(id));
	end
	tmp = reale(tmp);
	st.sn = st.sn(:) * tmp';
end
if length(Nd) > 1
	st.sn = reshape(st.sn, Nd);	% [(Nd)]
else
	st.sn = st.sn(:);	% [(Nd)]
end

%
% [J?,M] interpolation coefficient vectors.  will need kron of these later
%
for id=1:dd
	N = Nd(id);
	J = Jd(id);
	K = Kd(id);
	if isfield(st, 'kernel')
		[c, arg] = ...
		nufft_coef(om(:,id), J, K, st.kernel{id});	% [J?,M]
	else
		alpha = st.alpha{id};
		beta = st.beta{id};
		T = nufft_T(N, J, K, st.tol, alpha, beta);	% [J?,J?]
		[r, arg] = ...
		nufft_r(om(:,id), N, J, K, alpha, beta);	% [J?,M]
		c = T * r;
		clear T r
	end

	%
	% indices into oversampled FFT components
	%
	koff = nufft_offset(om(:,id), J, K);	% [M,1] to leftmost near nbr
	k0 = outer_sum([1:J]', koff');	% [J?,M] arbitrary integers
	kd{id} = mod(k0, K);	% [J?,M] {0,...,K?-1} (DFT indices)

	gam = 2*pi/K;

	switch st.phasing
	case {'real', 'none'}
		phase = 1;
	case 'complex'
		phase_scale = 1i * gam * (N-1)/2;
		phase = exp(phase_scale * arg);	% [J?,M] linear phase
	case 'flipreal'
		isodd = @(n) mod(n,2) == 1;
		phase = ones(size(k0)); % [J?,M]
		flip = isodd((kd{id} - k0) / K * (N-1)); % sign flip every K
		phase(flip) = -1; % sign flip (if N is even)
	otherwise
		fail('unknown phasing %s', st.phasing)
	end
	ud{id} = phase .* c;		% [J?,M]

end % id
clear c arg gam phase phase_scale koff k0 N J K

% st.M = size(om,1);
M = size(om,1);

%
% build sparse matrix that is [M,*Kd]
% with *Jd nonzero entries per frequency point
%
if st.printmem
	printm('Needs at least %g Gbyte RAM', prod(Jd)*M*8/2^30*2)
end

kk = kd{1};	% [J1,M]
uu = ud{1};	% [J1,M]
for id = 2:dd
	Jprod = prod(Jd(1:id));
	tmp = kd{id} * prod(Kd(1:(id-1)));
	kk = block_outer_sum(kk, tmp);		% outer sum of indices
	kk = reshape(kk, Jprod, M);
	uu = block_outer_prod(uu, ud{id});	% outer product of coefficients
	uu = reshape(uu, Jprod, M);
end % now kk and uu are [*Jd, M]

%
% handle phase shifts
%
uu = conj(uu); % [*Jd,M] ala Hermitian transpose of interpolation coefficients
switch st.phasing
case 'complex'
	phase = exp(1i * (om * st.n_shift(:))).';	% [1,M]
	uu = uu .* phase(ones(1,prod(Jd)),:);	% [*Jd,M]
	if streq(st.mode, 'table', 5) % moved from newfft_table_init.m to here
		st.phase_after = @(om) exp(1i * (om * col(st.n_shift))); % [M,1]
	end
case {'real', 'flipreal'}
	st.phase_before = newfft_phase_before(Kd, st.Nmid);
	if ~isempty(st.om) % precompute phase that goes after interpolation
		st.phase_after = newfft_phase_after(st.om, st.Nmid, st.n_shift);
	else
		st.phase_after = @(om) newfft_phase_after(om, st.Nmid, st.n_shift);
	end
case 'none'
	% do nothing
otherwise
	error 'bug'
end

mm = repmat(1:M, prod(Jd), 1); % [*Jd,M]
st.p = sparse(mm(:), 1+kk(:), uu(:), M, prod(Kd)); % [M, *Kd] sparse matrix
% sparse object, to better handle single precision operations!
st.p = Gsparse(st.p, 'odim', [M 1], 'idim', [prod(Kd) 1]);

st = strum(st, { ...
	'fft', @newfft_approx_for, '(x, [om])';
	'adj', @newfft_approx_adj, '(X, [om])';
	});


%
% in
%	x1	[J1,M]
%	x2	[J2,M]
% out
%	y	[J1,J2,M]	y(i1,i2,m) = x1(i1,m) + x2(i2,m)
%
function y = block_outer_sum(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 + xx2;			% [J1,J2,M]

function y = block_outer_prod(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 .* xx2;			% [J1,J2,M]


%
% newfft_phase_before()
% phase factor that gets multiplied by DFT (before interpolation)
%
function phase = newfft_phase_before(Kd, Nmid)

phase = 0;
for id = 1:length(Kd)
	tmp = 2 * pi * [0:Kd(id)-1] / Kd(id) * Nmid(id);
	phase = outer_sum(phase, tmp); % [(Kd)] when done
end
phase = exp(1i * phase);


%
% newfft_phase_after()
% phase factor that multiplies the DTFT (after interpolation)
%
function phase = newfft_phase_after(om, Nmid, n_shift)

phase = exp(1i * (om * col(n_shift - Nmid))); % [M,1]


%
% newfft_scale_tri()
% scale factors when kernel is 'linear'
% tri(u/J) <-> J sinc^2(J x)
%
function sn = newfft_scale_tri(N, J, K, Nmid)
nc = [0:N-1] - Nmid;
fun = @(x) J * nufft_sinc(J * x / K).^2;
cent = fun(nc);
sn = 1 ./ cent;

% try the optimal formula
tmp = 0;
LL = 3;
for ll=-LL:LL
	tmp = tmp + abs(fun(nc - ll*K)).^2;
end
sn = cent ./ tmp;


%
% newfft_test_time()
% compare compute times of real vs complex, sparse vs table
%
function newfft_test_time

Nd = [1 1] * 2^8;
[tmp om] = mri_trajectory('radial', {}, Nd, Nd);
pr length(om)
rng(0)
x0 = rand([Nd 1]);

arg = {om, Nd, 'dotime', ' '};
s0r = newfft(arg{:}, 'mode', 'table0');
s0c = newfft(arg{:}, 'mode', 'table0', 'phasing', 'complex');
s1r = newfft(arg{:}, 'mode', 'table1');
ssr = newfft(arg{:}, 'mode', 'sparse');
ssc = newfft(arg{:}, 'mode', 'sparse', 'phasing', 'complex');

tmp = @(st) newfft_test_time_one(st, x0, 2);
tmp(s0r)
tmp(s0c)
tmp(s1r)
tmp(ssr)
tmp(ssc)

function newfft_test_time_one(st, x0, nn)
st.fft(x0); % warm up
tmp = [st.mode ' ' st.phasing(1)];
cpu etic
for ii=1:nn
	st.fft(x0); % trial
end
cpu('etoc', tmp)


%
% newfft_test()
%
function newfft_test

newfft_test_time
prompt

modes = {'sparse', 'table0', 'table1'};

phasings = {'complex', 'real'};
% 'flipreal' no longer needed thanks to floor(N/2)
% 'none' is for internal table only

Nd_list = [20 10 8];
for id=0:3
	if id == 0
		Nd = 1 + Nd_list(1); % test odd case
	else
		Nd = Nd_list(1:id);
	end
	dd = length(Nd);

	rng(0)
	om = 3 * 2 * pi * (rand(100,length(Nd)) - 0.5);
%	om = sort(om);
%	om = linspace(-1,1,601)' * 3*pi;
%	om = [-8:8]'/2 * pi;
	%om = 'epi';
	x0 = rand([Nd 1]);
%	x0 = zeros([Nd 1]); x0(1) = 1; % unit vector for testing
%	x0 = ones([Nd 1]);

	% exact
	st_e = newfft(om, Nd, 'mode', 'exact');
	Xe = st_e.fft(x0);
	xe = st_e.adj(Xe);
	if 0
		pe = st_e.p;
		equivs(pe * x0(:), Xe)
		equivs(reshape(pe' * Xe, [Nd 1]), xe)
	end

	ktypes = {{'linear', 'Jd', 2*ones(1,dd)}, ...
		'minmax:unif', ...
		{'minmax:user', 'alpha', num2cell(ones(1,dd)), ...
			'beta', num2cell(0.5 * ones(1,dd))}, ...
		{'diric', 'Jd', 2*Nd-0, 'oversample', []}, ...
		'minmax:kb', 'minmax:tuned', ...
		'kb:minmax', 'kb:beatty', ...
		{'kb:user', 'kb_m', 0*Nd, 'kb_alf', 2.34 * 6 + 0*Nd}
		};

	for ii=4:length(ktypes) % skip poor ones
		ktype = ktypes{ii};
		if ~iscell(ktype), ktype = {ktype}; end

	for jj=1:length(modes)
	for ip=1:length(phasings)
		sc = newfft(st_e.om, st_e.Nd, 'phasing', 'complex', ...
			'mode', 'table0', 'ktype', ktype{:});

		st = newfft(st_e.om, st_e.Nd, 'phasing', phasings{ip}, ...
			'mode', modes{jj}, 'ktype', ktype{:});

		if streq(st.phasing, 'complex') && streq(st.mode, 'table1')
			continue
		end

%		pr minmax(st.sn)
		pad = @(s,n) [s blanks(n-length(s))];
		key = [sprintf('%2d ', st.Jd(1)) st.ktype];
		key = [st.mode ' ' num2str(id) st.phasing(1) ' ' pad(key,16)];
		Xs = st.fft(x0);
		max_percent_diff(Xe, Xs, key)
%		plot(abs(Xe), abs(Xs), 'o'), prompt

		xs = st.adj(Xe);
		max_percent_diff(xe, xs, key)
	end % ip
	end % jj
	end % ii
end
