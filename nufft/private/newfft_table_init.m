 function st = newfft_table_init(st, varargin)
%function st = newfft_table_init(st, varargin)
%|
%| Initialize structure for d-dimension NUFFT using table-based interpolator,
%| This should be called only by newfft for its 'table0' or 'table1' mode!
%| Note: default oversample factor is 2^11 or 2^13
%|
%| in
%|	st		structure initialized in newfft.m
%|			uses st.phasing to choose complex or real kernel
%|
%| option
%|	'how'	char	'slow' | 'ratio' | 'fast'
%|			(default depends on kernel type)
%|
%| out
%|	st	strum	(see newfft)
%|
%| Copyright 2008-7-11, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

arg.how = '';
arg = vararg_pair(arg, varargin);

if isempty(arg.how)
	if streq(st.ktype, 'kb:', 3) % todo: others?
		arg.how = 'ratio';
	else
		arg.how = 'fast';
	end
end

st.phase_shift = []; % compute on-the-fly

Nd = st.Nd;
Jd = st.Jd;
Kd = st.Kd;
dd = length(Nd);
Ld = st.oversample;
if isempty(Ld)
	switch st.mode
	case 'table0'
		st.order = 0; % 0th order
		Ld = 2^11; % default
	case 'table1'
		st.order = 1; % 1st order
		Ld = 2^9; % default
	otherwise
		fail 'bad mode'
	end
end

% dimensionality of input space (usually 2 or 3)
if dd > 1 && length(Ld) == 1
	Ld = repmat(Ld, [1 dd]); % allow scalar L to be used for all dimensions
end
st.oversample = Ld;
st.Ld = Ld;

if dd ~= length(Jd) || dd ~= length(Kd) || dd ~= length(Ld)
	whos
	fail 'inconsistent dim'
end

% desired scaling factors, using fake omega=[0]
ktype_args = {'ktype', st.ktype}; % user's interpolator
if ~isempty(st.kb_m)
	ktype_args = {ktype_args{:}, 'kb_m', st.kb_m, 'kb_alf', st.kb_alf};
elseif ~isempty(st.alpha)
	ktype_args = {ktype_args{:}, 'alpha', st.alpha, 'beta', st.beta};
end

tmp = newfft(zeros(size(Nd)), st.Nd, 'Jd', st.Jd, 'Kd', st.Kd, ...
	'n_shift', 0*st.n_shift, ktype_args{:}, 'mode', 'sparse', ...
	'phasing', st.phasing);
st.sn = tmp.sn;

if streq(st.phasing, 'flipreal')
	iseven = @(n) mod(n,2) == 0;
	st.flips = iseven(st.Nd); % [1 dd]
end


% 1D tables for each dimension
for id = 1:dd

	ktype_args = {'ktype', st.ktype}; ... % user's interpolator
	if ~isempty(st.kb_m)
		ktype_args = {ktype_args{:}, ...
			'kb_m', st.kb_m(id), 'kb_alf', st.kb_alf(id)};
	elseif ~isempty(st.alpha)
		ktype_args = {ktype_args{:}, ...
			'alpha', {st.alpha{id}}, 'beta', {st.beta{id}}};
	end

	make_args = {Nd(id), Jd(id), Kd(id), Ld(id), ktype_args, st.phasing};
	st.h{id} = newfft_table_make(arg.how, make_args{:});

	if 0 % test fast vs slow
		[h0 t0] = newfft_table_make('slow', make_args{:});
		h = st.h{id};
		jf plc 2 2
		jf sub 1, plot(t0, real(h0), 'c-', t0, real(h), 'y.')
		jf sub 2, plot(t0, real(h0-h), 'y.')
		jf sub 3, plot(t0, imag(h0), 'c-', t0, imag(h), 'y.')
		jf sub 4, plot(t0, imag(h0-h), 'y.')
		max_percent_diff(h0, h)
	prompt
	end

	switch st.phasing
	case 'complex'
		if isreal(st.h{id})
			warn 'real kernel?'
%			st.h{id} = complex(st.h{id});
		end

	case {'real', 'flipreal', 'none'}
		if ~isreal(st.h{id})
			fail 'nonreal kernel bug'
		end

	otherwise
		fail bug
	end

end

st = strum(st, { ...
	'fft', @newfft_approx_for, '(x, [om])';
	'adj', @newfft_approx_adj, '(X, [om])';
	'interp_for', @newfft_table_for, '(Xk, [om])';
	'interp_adj', @newfft_table_adj, '(Xk, [om])';
	'plot', @newfft_table_plot, '()';
	});


% newfft_table_make()
function [h, t0] = newfft_table_make(how, N, J, K, L, ktype_args, phasing)

if streq(phasing, 'flipreal')
	phasing = 'none'; % trick: get pure real kernel (with no sign flips)
end
% note: if phasing is already 'real' then keep it that way!

mode_args = {'mode', 'sparse', 'Jd', J, 'n_shift', 0, 'phasing', phasing};

t0 = [-J*L/2:J*L/2]' / L; % [J*L+1]

switch how

% This is a slow and inefficient (but simple) way to get the table
% because it builds a huge sparse matrix but only uses 1 column!
case 'slow'

	om0 = t0 * 2*pi/K; % gam
	s1 = newfft(om0, N, 'Kd', K, ktype_args{:}, mode_args{:});
	h = full(s1.p(:,1));

% This way is "J times faster" than the slow way, but still not ideal.
% It works for any user-specified interpolator.
case 'fast'

	t1 = J/2-1 + [0:(L-1)]' / L;	% [L]
	om1 = t1 * 2*pi/K;		% * gam
	s1 = newfft(om1, N, 'Kd', K, ktype_args{:}, mode_args{:});
	h = col(full(s1.p(:,J:-1:1))); % [J*L+0]
%	h = [h; 0]; % [J*L+1]
	h = [h; h(1)]; % [J*L+1] - assuming symmetric

% This efficient way uses only "J" columns of sparse matrix!
% The trick to this is to use fake small values for N and K,
% which works for interpolators that depend only on the ratio K/N.
case 'ratio' % e.g., 'minmax:kb' | 'kb:*'

	Nfake = J;
	Kfake = Nfake * K/N;
	t1 = J/2-1 + [0:(L-1)]' / L;	% [L]
	om1 = t1 * 2*pi/Kfake;		% "gam"
	s1 = newfft(om1, Nfake, 'Kd', Kfake, ktype_args{:}, mode_args{:});
	h = col(full(s1.p(:,J:-1:1))); % [J*L+0]
%	h = [h; 0]; % [J*L+1]
	h = [h; h(1)]; % [J*L+1] assuming symmetric

	if streq(phasing, 'complex')
		h = h .* exp(1i * pi * t0 * (1/K - 1/Kfake)); % fix phase
	end

otherwise
	fail('bad type %s', how)
end


% newfft_table_plot()
function newfft_table_plot(st)
color = 'cym';
arg = {};
for id = 1:st.dd
	J = st.Jd(id);
	L = st.Ld(id);
	t0 = [-J*L/2:J*L/2]' / L; % [J*L+1]
	arg = {arg{:}, t0, real(st.h{id}), [color(id) '-']};
	arg = {arg{:}, t0, imag(st.h{id}), [color(id) '--']};
end
plot(arg{:})


% newfft_table_for()
function X = newfft_table_for(st, Xk, om)
if ~isvar('om') || isempty(om)
	om = st.om;
end
if isempty(om)
	fail 'need om or st.om'
end
if st.dd ~= size(om,2), fail('omega needs %d columns', dd), end

X = nufft_table_interp(st, Xk, st.order, st.flips, om);


% newfft_exact_adj()
function Xk = newfft_table_adj(st, X, om)
if ~isvar('om') || isempty(om)
	om = st.om;
end
if isempty(om)
	fail 'need om or st.om'
end
if st.dd ~= size(om,2), fail('omega needs %d columns', dd), end

Xk = nufft_table_adj(st, X, st.order, st.flips, om);
