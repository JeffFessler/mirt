 function st = nufft_init(om, Nd, Jd, Kd, varargin)
%function st = nufft_init(om, Nd, Jd, Kd, [n_shift,] ...)
%|
%| Initialize structure for d-dimension NUFFT using KB interpolator,
%| particularly the interpolation matrix in sparse format.
%| caution: this routine can require a lot of memory!
%| in
%|	om [M d]	"digital" frequencies in radians
%|	Nd [d]		image dimensions (N1,N2,...,Nd)
%|	Jd [d]		# of neighbors used (in each direction)
%|	Kd [d]		FFT sizes (should be >= N1,N2,...)
%|
%| options
%|	n_shift [d]	n = 0-n_shift to N-1-n_shift (must be first)
%|	'minmax:kb'	minmax interpolator with excellent KB scaling!
%|				(minmax:kb is recommended, and used by default)
%|	'minmax:tuned'	minmax interpolator, somewhat numerically tuned
%|	'minmax:user'	minmax interpolator with user ({alpha}, {beta})
%|	'uniform'	uniform scaling factors (not recommended)
%|	'kaiser'	kaiser-bessel (KB) interpolator (minmax best alpha, m)
%|			or 'kaiser', [alpha], [m] to specify parameter (vectors)
%|	'linear'	linear interpolator (a terrible straw man)
%|	kernel		user-provided function handle interpolation kernel(k,J)
%|			(or a cell array of kernels, one for each dimension)
%|	'table'		use table-based interpolation rather than sparse matrix.
%|			this can save a lot of memory for large problems.
%|			example ..., 'table', 2^11, 'minmax:kb'
%|			where 2^11 is the table over-sampling factor.
%|
%| out
%|	st.p		[M *Kd]		sparse interpolation matrix
%|					(or empty if table-based)
%|	st.sn		[(Nd)]		scaling factors
%|	st.Nd,Jd,Kd,om	copies of inputs
%|
%| *Nd is shorthand for prod(Nd).
%| (Nd) is shorthand for (N1,N2,...,Nd)
%|
%| Like fft(), the NUFFT expects the signals to be x(0,0), ...
%| Use n_shift = [N1/2, N2/2, ...] for x(-N1/2,-N2/2,...), ...
%|
%| Copyright 2002-5-30, Jeff Fessler, University of Michigan

if nargin == 1 && streq(om, 'test'), nufft_init_test, return, end
if nargin < 4, ir_usage, end

% dimensionality of input space (usually 2 or 3)
dd = length(Nd);
if dd ~= length(Jd) || dd ~= length(Kd)
	fail 'inconsistent dim'
end
if any(Kd < Nd), warn 'Kd < Nd unlikely to work.  Try Kd=2*Nd', end
if any(Jd > Kd)
	Jd_new = min(Jd, Kd);
	fun = @(x) regexprep(num2str(x), '\s+', ' ');
	warn('Jd = [%s] > Kd = [%s]; changing to Jd = [%s]', ...
		fun(Jd), fun(Kd), fun(Jd_new))
	Jd = Jd_new; clear Jd_new fun
end

% special cases of input sampling pattern
if ischar(om)
	om = nufft_samples(om, Nd);
end
if dd ~= size(om,2), fail('omega needs %d columns', dd), end

%
% process optional arguments
%

% n_shift argument? (must be first)
if length(varargin) > 0 && isnumeric(varargin{1})
	n_shift = varargin{1};
	if dd ~= length(n_shift)
		fail('n_shift needs %d columns', dd)
	end
	varargin = {varargin{2:end}};
else
	n_shift = zeros(size(Nd));
end
st.n_shift = n_shift;

% default/recommended interpolator is minmax with KB scaling factors
if length(varargin) == 0
	varargin = {'minmax:kb'};
end

st.alpha = {};
st.beta = {};
is_kaiser_scale = false;

% table based?
if ischar(varargin{1}) && streq(varargin{1}, 'table')
	st = nufft_table_init(om, Nd, Jd, Kd, n_shift, varargin{2:end});
	return
end

ktype = varargin{1};

% cell array of kernel functions: {kernel1, kernel2, ..., kernelD}
if isa(ktype, 'cell')
	if isa(ktype{1}, 'inline') || isa(ktype{1}, 'function_handle')
		ktype = 'function_handle';
		if length(varargin) > 1, fail 'excess arguments?', end
		if length(varargin{1}) ~= dd, fail 'wrong # of kernels', end
		st.kernel = varargin{1};
	else
		fail 'cell array should be function_handles!?'
	end

% or a single function handle for all dimension
elseif isa(ktype, 'inline') || isa(ktype, 'function_handle')
	ktype = 'function_handle';
	if length(varargin) > 1, fail 'excess arguments?', end
	for id = 1:dd
		st.kernel{id} = varargin{1}; % all same
	end

% or a string that describes the type of interpolator
elseif ~ischar(ktype)
	fail 'non-string kernel type?'

end

st.ktype = ktype;

% set up whatever is needed for each interpolator
switch ktype
case 'function_handle'
	% already did it above

% linear interpolator straw man
case 'linear'
	ktype = 'function_handle';
	kernel = @(k,J) (1 - abs(k/(J/2))) .* (abs(k) < J/2);
	for id = 1:dd
		st.kernel{id} = kernel;
	end

% KB interpolator
case 'kaiser'
	is_kaiser_scale = true;

	% with minmax-optimized parameters
	if length(varargin) == 1
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id));
		end

	% with user-defined parameters
	elseif length(varargin) == 3
		alpha_list = varargin{2};
		m_list = varargin{3};
		if (length(alpha_list) ~= dd) || (length(m_list) ~= dd)
			fail('#alpha=%d #m=%d vs dd=%d', ...
				length(alpha_list), length(m_list), dd)
		end
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id), ...
					alpha_list(id), m_list(id));
		end
	else
		fail 'kaiser should have no arguments, or both alpha and m'
	end

% minmax interpolator with KB scaling factors (recommended default)
case 'minmax:kb'
	for id = 1:dd
		[st.alpha{id}, st.beta{id}] = ...
			nufft_alpha_kb_fit(Nd(id), Jd(id), Kd(id));
	end

% minmax interpolator with numerically "tuned" scaling factors
case 'minmax:tuned'
	for id = 1:dd
		[st.alpha{id}, st.beta{id}, ok] = ...
			nufft_best_alpha(Jd(id), 0, Kd(id)/Nd(id));
		if ~ok, fail 'unknown J,K/N', end
	end

% minmax interpolator with user-provided scaling factors
case 'minmax:user'
	if length(varargin) ~= 3, fail 'user must provide alpha/beta', end
	st.alpha = varargin{2};
	st.beta = varargin{3};
	if length(st.alpha) ~= dd || length(st.beta) ~= dd
		fail 'alpha/beta size mismatch'
	end

case 'uniform'
	for id = 1:dd
		st.alpha{id} = 1;
		st.beta{id} = 0;
	end

otherwise
	fail('unknown kernel type %s', ktype)
end

nufft_check_dft(om, Nd, ktype)

st.tol = 0;

st.Jd = Jd;
st.Nd = Nd;
st.Kd = Kd;

M = size(om,1);
st.M = M;
st.om = om;

% scaling factors: "outer product" of 1D vectors
st.sn = 1;
for id=1:dd
	if is_kaiser_scale
		nc = [0:Nd(id)-1]'-(Nd(id)-1)/2;
		tmp = 1 ./ kaiser_bessel_ft(...
			nc/Kd(id), Jd(id), st.kb_alf(id), st.kb_m(id), 1);
	elseif streq(ktype, 'function_handle')
		tmp = 1 ./ nufft_interp_zn(0, Nd(id), Jd(id), Kd(id), st.kernel{id});
	else
		tmp = nufft_scale(Nd(id), Kd(id), st.alpha{id}, st.beta{id});
	end
	st.sn = st.sn(:) * tmp';
end
if length(Nd) > 1
	st.sn = reshape(st.sn, Nd); % [(Nd)]
else
	st.sn = st.sn(:); % [(Nd)]
end

% [J? M] interpolation coefficient vectors.  will need kron of these later
for id=1:dd
	N = Nd(id);
	J = Jd(id);
	K = Kd(id);
	if isfield(st, 'kernel')
		[c, arg] = ...
		nufft_coef(om(:,id), J, K, st.kernel{id}); % [J? M]
	else
		alpha = st.alpha{id};
		beta = st.beta{id};
		T = nufft_T(N, J, K, st.tol, alpha, beta); % [J? J?]
		[r, arg] = ...
		nufft_r(om(:,id), N, J, K, alpha, beta); % [J? M]
		c = T * r;	clear T r
	end

	gam = 2*pi/K;
	phase_scale = 1i * gam * (N-1)/2;

	phase = exp(phase_scale * arg); % [J? M] linear phase
	ud{id} = phase .* c; % [J? M]

	% indices into oversampled FFT components
	koff = nufft_offset(om(:,id), J, K); % [M 1] to leftmost near nbr
	kd{id} = mod(outer_sum([1:J]', koff'), K) + 1; % [J? M] {1,...,K?}
	if id > 1 % trick: pre-convert these indices into offsets!
		kd{id} = (kd{id}-1) * prod(Kd(1:(id-1)));
	end

end, clear c arg gam phase phase_scale koff N J K

% build sparse matrix that is [M *Kd]
% with *Jd nonzero entries per frequency point
if dd >= 3
	tmp = prod(Jd)*M*8/10^9*2;
	if tmp > 1. % only display if more than 1GB
		printm('Needs at least %g Gbyte RAM', tmp)
	end
end

kk = kd{1}; % [J1 M]
uu = ud{1}; % [J1 M]
for id = 2:dd
	Jprod = prod(Jd(1:id));
	kk = block_outer_sum(kk, kd{id}); % outer sum of indices
	kk = reshape(kk, Jprod, M);
	uu = block_outer_prod(uu, ud{id}); % outer product of coefficients
	uu = reshape(uu, Jprod, M);
end % now kk and uu are [*Jd M]

% apply phase shift
% pre-do Hermitian transpose of interpolation coefficients
phase = exp(1i * (om * n_shift(:))).'; % [1 M]
uu = conj(uu) .* phase(ones(1,prod(Jd)),:); % [*Jd M]

mm = [1:M]; mm = mm(ones(prod(Jd),1),:); % [*Jd M]
% make sparse matrix, ensuring arguments are double for stupid matlab
st.p = sparse(mm(:), double(kk(:)), double(uu(:)), M, prod(Kd));
% sparse object, to better handle single precision operations!
st.p = Gsparse(st.p, 'odim', [M 1], 'idim', [prod(Kd) 1]);


% block_outer_sum()
%
% in
%	x1	[J1 M]
%	x2	[J2 M]
% out
%	y	[J1 J2 M]	y(i1,i2,m) = x1(i1,m) + x2(i2,m)
%
function y = block_outer_sum(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]); % [J1 1 M] from [J1 M]
xx1 = xx1(:,ones(J2,1),:); % [J1 J2 M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]); % [1 J2 M] from [J2 M]
xx2 = xx2(ones(J1,1),:,:); % [J1 J2 M], emulating ndgrid
y = xx1 + xx2; % [J1 J2 M]


% block_outer_prod()
function y = block_outer_prod(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]); % [J1 1 M] from [J1 M]
xx1 = xx1(:,ones(J2,1),:); % [J1 J2 M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]); % [1 J2 M] from [J2 M]
xx2 = xx2(ones(J1,1),:,:); % [J1 J2 M], emulating ndgrid
y = xx1 .* xx2; % [J1 J2 M]


% nufft_check_dft()
% see if they are DFT samples; if so, give warning
function nufft_check_dft(om, Nd, ktype)
kk = om / (2*pi) .* repmat(Nd(:)', [nrow(om) 1]);
tol = 1e-6;
if all(col(abs(round(kk) - kk)) < tol) ...
	&& any(om(:)) ... % trick: ignore all-zero om used in nufft_table_init()
	&& (streq(ktype, 'minmax:kb') || streq(ktype, 'kaiser'))
	warn('DFT samples with ktype="%s" has suboptimal accuracy', ktype)
end


% nufft_init_test()
% for a more complete test, use "nufft test"
function nufft_init_test
Nd = [20 10];
st = nufft_init('epi', Nd, [5 5], 2*Nd);
if 0
	clf
	om = st.om;
	plot(om(:,1), om(:,2), 'o-')
	axis_pipi set
end
st;
st.alpha{1};
