 function X = nufft(x, st)
%function X = nufft(x, st)
%|
%| Compute d-dimensional NUFFT of signal/image x
%|
%| in
%|	x	[N1 N2 ... Nd (L)]	L input image(s) of size
%|						N1 x N2 x ... x Nd
%|	st	structure		precomputed by nufft_init()
%| out
%|	X	[M (L)]			output spectra
%|
%| Copyright 2003-5-30, Jeff Fessler, University of Michigan

% if no arguments, then run some self tests
if (nargin < 1) || (nargin == 1 && streq(x, 'test'))
	help([mfilename '.m'])

	printm('Starting nufft() self test; after a wait, shows small errors')
	if 0
		Jd = [5 4 4];
		Nd = [23 15 19];
		alpha_user = {1, 1, 1}; % default alpha, beta
		beta_user = {0.5, 0.5, 0.5};
	else
		Jd = [5 6];
		Nd = [60 75];
		alpha_user = {1, 1};
		beta_user = {0.5, 0.5};
	end
	Kd = 2 * Nd;
	gam = 2*pi ./ Kd;
	n_shift = zeros(size(Nd));

	if 1
		oldarg = {Nd(1), Nd(2), Jd(1), Jd(2), Kd(1), Kd(2)};
		printm('err alf1 %g best %g', ...
			max(col(nufft2_err_mm('all', oldarg{:}, 1))), ...
			max(col(nufft2_err_mm('all', oldarg{:}, 'best'))) )
	end
	rng(0)
	x = randn(Nd);
%	x = [[1:Nd(1)]'*ones(1,3), ones(Nd(1), Nd(2)-3)]; % test signal
%	x = zeros(N1, N2); x(1,1) = 1;
	if 0 % test with uniform frequency locations
		o1 = 2 * pi * [0:(N1-1)]' / N1;
		o2 = 2 * pi * [0:(N2-1)]' / N2;
		[o1, o2] = ndgrid(o1, o2);
		Xf = fft2(x);
	else
		if length(Nd) == 3 % nonuniform frequencies
			[o1, o2, o3] = ndgrid(linspace(0,gam(1),11), ...
				linspace(0,gam(2),13), linspace(0,gam(3),5));
			om = [	[o1(:); [0 7.2 2.6 3.3]'], ...
				[o2(:); [0 4.2 -1. 5.5]'], ...
				[o3(:); [0 1.1 -2. 3.4]'] ];
		else
			[o1, o2] = ndgrid(linspace(-3*gam(1),3*gam(1),41), ...
				linspace(-2*gam(2),gam(2),31));
			om = [	[o1(:); [0 7.2 2.6 3.3]'], ...
				[o2(:); [0 4.2 -1. 5.5]'] ];
		end
		Y.d = dtft(x, om, 'n_shift', n_shift);
	end

	s.tab = nufft_init(om, Nd, Jd, Kd, n_shift, 'table', 2^12, 'minmax:kb');
	Y.tab = nufft(x, s.tab);
	printm('table0		max%%diff = %g', max_percent_diff(Y.d, Y.tab))

	s.mmkb = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:kb');
	Y.mmkb = nufft(x, s.mmkb);
	printm('minmax:kb	max%%diff = %g', max_percent_diff(Y.d, Y.mmkb))

	if 0 % test multiple input case
		xx = x; xx(:,:,:,2) = x; xx(:,:,:,3) = x;
		Y.tmp = nufft(xx, s.mmkb);
		Y.tmp = Y.tmp(:,1);
		printm('multi	max%%diff = %g', max_percent_diff(Y.mmkb, Y.tmp))
	return
	end

	% kaiser with minmax best alpha,m
	s.kb = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser');
	Y.kb = nufft(x, s.kb);
	printm('kaiser		max%%diff = %g', max_percent_diff(Y.d, Y.kb))

	% kaiser with user-specified supoptimal alpha,m for comparison
	s.ku = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser', ...
		s.kb.kb_alf + 0.1*ones(size(s.kb.kb_alf)), s.kb.kb_m);
	Y.ku = nufft(x, s.ku);
	printm('kaiser-user	max%%diff = %g', max_percent_diff(Y.d, Y.ku))

	s.mmtu = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:tuned');
	Y.mmtu = nufft(x, s.mmtu);
	printm('minmax:tuned	max%%diff = %g', max_percent_diff(Y.d, Y.mmtu))

	s.mm = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:user', ...
		alpha_user, beta_user);
	Y.mm = nufft(x, s.mm);
	printm('minmax:user	max%%diff = %g', max_percent_diff(Y.d, Y.mm))

	if 0 % test 'uniform' scaling
		s.un = nufft_init(om, Nd, Jd, Kd, n_shift, 'uniform')
		Y.un = nufft(x, s.un);
		printm('user-unif max%%diff = %g', max_percent_diff(Y.mm, Y.un))
	return
	end

%	s1 = nufft_init(om, Nd, Jd, Kd, n_shift, 'loop', 'kaiser');
%	printm('kb loop max %% difference = %g', max_percent_diff(sn.p,s1.p))

 if 0 % test against old 2D version
	s2 = nufft2_init_kb(om, oldarg{:}, n_shift, 'kaiser');
%	Y2 = nufft2(x, s2);
%	printm('KB 2d vs nd max%%diff = %g', max_percent_diff(Y2, Y.kb))
	printm('KB 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.kb.p))

	s2 = nufft2_init(om, oldarg{:}, n_shift, 0, 'best');
	printm('tuned 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.mmtu.p))

	s2 = nufft2_init(om, oldarg{:}, n_shift, 0); % minmax
	printm('user 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.mm.p))

 end
return
end

Nd = st.Nd;
Kd = st.Kd;

dims = size(x);
dd = length(Nd);
if ndims(x) < dd, fail 'input signal has too few dimensions', end
if any(dims(1:dd) ~= Nd), fail 'input signal has wrong size', end


% the usual case is where L=1, i.e., there is just one input signal.
if ndims(x) == dd
	x = x .* st.sn; % apply scaling factors
	Xk = col(fftn_fast(x, Kd)); % [*Kd] oversampled FFT, padded at end

% otherwise, collapse all excess dimensions into just one
else
	xx = reshape(x, [prod(Nd) prod(dims((dd+1):end))]); % [*Nd *L]
	L = size(xx, 2);
	Xk = zeros(prod(Kd), L, class(xx)); % [*Kd *L]
	for ll=1:L
		xl = reshape(xx(:,ll), [Nd 1]); % l'th signal
		xl = xl .* st.sn; % scaling factors
		Xk(:,ll) = col(fftn_fast(xl, [Kd 1]));
	end
end


% interpolate using precomputed sparse matrix
% or with tabulated interpolator
if ~isfield(st, 'interp_table')
	X = st.p * Xk; % [M *L]
else
	X = feval(st.interp_table, st, Xk);
end

if ndims(x) > dd
	X = reshape(X, [st.M dims((dd+1):end)]); % [M (L)]
end
