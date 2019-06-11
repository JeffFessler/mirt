 function st = nufft_table_init(om, Nd, Jd, Kd, n_shift, Ld, varargin)
%function st = nufft_table_init(om, Nd, Jd, Kd, n_shift, Ld, ...)
%|
%| Initialize structure for d-dimension NUFFT using table-based interpolator,
%| This should be called only by nufft_init with its 'table' argument!
%| (So this m-file belongs in the private directory.)
%|
%| in
%|	om, Nd, Jd, Kd, n_shift		see nufft_init.m
%|	Ld [d]				table over-sampling factor(s)
%|
%| optional arguments (same as nufft_init.m, to specify interp type)
%| out
%|	st.interp_table		function: X = st.interp_table(st, Xk);
%|	st.interp_table_adj	function: Xk = st.interp_table_adj(st, Xk);
%|	st.phase_shift		phase shift due to n_shift
%|	st.*			same as in nufft_init.m
%|
%| Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

if nargin < 5, help(mfilename), error args, end

% dimensionality of input space (usually 2 or 3)
dd = length(Nd);
if dd > 1 && length(Ld) == 1
	Ld = Ld * ones(1,dd); % allow scalar L to be used for all dimensions
end
if dd ~= length(Jd) || dd ~= length(Kd) || dd ~= length(Ld)
	printm('dd:'), disp(dd)
	printm('Jd:'), disp(Jd)
	printm('Kd:'), disp(Kd)
	printm('Ld:'), disp(Ld)
	error 'inconsistent dim'
end
if dd ~= size(om,2), error('omega needs %d columns', dd), end

% D-dim NUFFT structure with desired scaling factors, using fake omega=0
st = nufft_init(zeros(size(Nd)), Nd, Jd, Kd, 0*n_shift, varargin{:});
st.p = [];
st.om = om;
st.Ld = Ld;
st.M = size(om,1);
st.n_shift = n_shift; % trick: store it in structure for history

% phase shift associated with indexing n_shift
% trick: because we put this phase shift here, we use "0" for n_shift below
% when building the tables.  in fact, n_shift definitely cannot be in table.
if any(n_shift ~= 0)
	st.phase_shift = exp(1i * (om * n_shift(:))); % [M 1]
end

% 1D tables for each dimension
for id = 1:dd
	J = Jd(id);
	L = Ld(id);
	K = Kd(id);

	% This is a slow and inefficient way to get the table
	% since it builds a huge sparse matrix but only uses 1 column!
	t0 = [-J*L/2:J*L/2]' / L; % [J*L+1]
	if 0
		om0 = t0 * 2*pi/K; % gam(id)
		s1 = nufft_init(om0, Nd(id), J, K, ...
			0*n_shift(id), varargin{:}); % user's interpolator
		h0 = full(s1.p(:,1));
	end

	% This efficient way uses only "J" columns of sparse matrix!
	% The trick to this is to use fake small values for N and K,
	% which works for interpolators that depend only on the ratio K/N.
	ktype = varargin{1};
	if 0 && (streq(ktype, 'minmax:kb') || streq(ktype, 'kaiser'))

		if streq(ktype, 'minmax:kb')
			args = {'minmax:user', {st.alpha{id}}, {st.beta{id}}};
		elseif streq(ktype, 'kaiser')
			args = {'kaiser'};
			if length(varargin) > 1, warning 'bug?', end
		else
			error 'bug'
		end

		Nfake = J;
		Kfake = Nfake * K/Nd(id);
		t1 = J/2-1 + [0:(L-1)]' / L;	% [L]
		om1 = t1 * 2*pi/Kfake;		% "gam"
		s1 = nufft_init(om1, Nfake, J, Kfake, 0, args{:});

		h = [];
		for jj=J:-1:1
			h = [h; full(s1.p(:,jj))]; % stack up pieces
		end
		h = [h; 0]; % finally J*L+1 long.
		h = h .* exp(1i * pi * t0 * (1/K - 1/Kfake)); % fix phase

		if 0	% testing
			clf, subplot(211)
			plot(t0, real(h0), 'c-', t0, real(h), 'y.')
			plot(t0, real(h0-h), 'y.')
			subplot(212)
			plot(t0, imag(h0), 'c-', t0, imag(h), 'y.')
			plot(t0, imag(h0-h), 'y.')
			printf('diff=%g', max_percent_diff(h0, h))
			keyboard
		end
		st.h{id} = double(h); % table mex files want double


	% This way is "J times faster" than the slow way, but still not ideal.
	% It works for any user-specified interpolator.
	else
		t1 = J/2-1 + [0:(L-1)]' / L; % [L]
		om1 = t1 * 2*pi/K; % gam(id)

		% trick: handle {'kaiser', [kb_alf], [kb_m]} in multidim case
		if streq(varargin{1}, 'kaiser') && length(varargin) == 3
			kb_alf = varargin{2};
			kb_m = varargin{3};
			if any(kb_alf ~= kb_alf(1)) || any(kb_m ~= kb_m(1))
				error 'table based needs identical interpolator for each dimension'
			end
			varargin{2} = kb_alf(1);
			varargin{3} = kb_m(1);
		end

		s1 = nufft_init(om1, Nd(id), J, K, 0*n_shift(id), varargin{:});

		h = [];
		for jj=J:-1:1
			h = [h; full(s1.p(:,jj))];
		end
		h = [h; 0];

		st.h{id} = double(h); % table mex files want double
		if 0	% testing
			clf, subplot(211)
			plot(t0, real(h0), 'c-', t0, real(h), 'y.')
			subplot(212)
			plot(t0, imag(h0), 'c-', t0, imag(h), 'y.')
			printf('diff=%g', max_percent_diff(h0, h))
			keyboard
		end
	end

	if isreal(st.h{id})
		if J == 1
			st.h{id} = complex(st.h{id});
		else
			whos
			warning 'real kernel?'
			keyboard
		end
	end

end

% interface routines
% st.interp_table = inline('nufft_table_interp(st, Xk)', 'st', 'Xk');
st.interp_table		= @nufft_table_interp;
st.interp_table_adj	= @nufft_table_adj;
