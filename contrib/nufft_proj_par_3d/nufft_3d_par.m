 function X = nufft_3d_par(x, st, om_xy)
%function X = nufft_3d_par(x, st, om_xy)
% Compute 1D NUFFT + 2D NUFFT of signal/image x
% 1D NUFFT along Wz direction
% 2D NUFFT on Wx-Wy plane
% in
%	x	[N1,N2,...,Nd,(3)]	3D input image(s) of size
%						N1 x N2 x ... x Nd
%	st	structure		precomputed by nufft_init()
%   st.xy               interpolation parameters for Wx, Wy coordinates 
%   st.z                interpolation parameters for Wz coordinates 
%   om_xy               frequency locations on Wx-Wy plane  
% out
%	X	[M,(3)]			output spectra
%
% Bsed on nufft
% Copyright 2003-5-30, Jeff Fessler, The University of Michigan
% 2007-6-7 modified for Fourier-based forward projection of 3D parallel 
% beam case(2D + 1D NUFFT), Yong Long

%
% if no arguments, then run some self tests
%
if (nargin < 1) || (nargin == 1 & streq(x, 'test'))
	help([mfilename '.m'])

	disp('Starting nufft() self test; after a wait, shows small errors')
	if 0
		Jd = [5 4 4];
		Nd = [23 15 19];
		alpha_user = {1, 1, 1};		% default alpha, beta
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
		printf('err alf1 %g best %g', ...
			max(col(nufft2_err_mm('all', oldarg{:}, 1))), ...
			max(col(nufft2_err_mm('all', oldarg{:}, 'best'))) )
	end
	rng(0)
	x = randn(Nd);
%	x = [[1:Nd(1)]'*ones(1,3), ones(Nd(1),Nd(2)-3)]; % test signal
%	x = zeros(N1,N2); x(1,1) = 1;
	if 0	% test with uniform frequency locations
		o1 = 2 * pi * [0:(N1-1)]' / N1;
		o2 = 2 * pi * [0:(N2-1)]' / N2;
		[o1, o2] = ndgrid(o1, o2);
		Xf = fft2(x);
	else
		if length(Nd) == 3	% nonuniform frequencies
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
		Y.d = dtft(x, om, n_shift);
	end

%	s.tab = nufft_init(om, Nd, Jd, Kd, 'table', 100, n_shift, 'minmax:kb');
%	Y.tab = nufft(x, s.tab);
%	printf('tab max%%diff = %g', max_percent_diff(Y.d, Y.tab))

	s.mmkb = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:kb');
	Y.mmkb = nufft(x, s.mmkb);
	printf('minmax:kb	max%%diff = %g', max_percent_diff(Y.d,Y.mmkb))

	if 0	% test multiple input case
		xx = x; xx(:,:,:,2) = x; xx(:,:,:,3) = x;
		Y.tmp = nufft(xx, s.mmkb);
		Y.tmp = Y.tmp(:,1);
		printf('multi	max%%diff = %g', max_percent_diff(Y.mmkb,Y.tmp))
	return
	end

	% kaiser with minmax best alpha,m
	s.kb = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser');
	Y.kb = nufft(x, s.kb);
	printf('kaiser		max%%diff = %g', max_percent_diff(Y.d,Y.kb))

	% kaiser with user-specified supoptimal alpha,m for comparison
	s.ku = nufft_init(om, Nd, Jd, Kd, n_shift, 'kaiser', ...
		s.kb.kb_alf + 0.1*ones(size(s.kb.kb_alf)), s.kb.kb_m);
	Y.ku = nufft(x, s.ku);
	printf('kaiser-user	max%%diff = %g', max_percent_diff(Y.d,Y.ku))

	s.mmtu = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:tuned');
	Y.mmtu = nufft(x, s.mmtu);
	printf('minmax:tuned	max%%diff = %g', max_percent_diff(Y.d,Y.mmtu))

	s.mm = nufft_init(om, Nd, Jd, Kd, n_shift, 'minmax:user', ...
		alpha_user, beta_user);
	Y.mm = nufft(x, s.mm);
	printf('minmax:user	max%%diff = %g', max_percent_diff(Y.d,Y.mm))

	if 0	% test 'uniform' scaling
		s.un = nufft_init(om, Nd, Jd, Kd, n_shift, 'uniform')
		Y.un = nufft(x, s.un);
		printf('user-unif max%%diff = %g', max_percent_diff(Y.mm,Y.un))
	return
	end

%	s1 = nufft_init(om, Nd, Jd, Kd, n_shift, 'loop', 'kaiser');
%	printf('kb loop max %% difference = %g', max_percent_diff(sn.p,s1.p))

 if 0	% test against old 2D version
	s2 = nufft2_init_kb(om, oldarg{:}, n_shift, 'kaiser');
%	Y2 = nufft2(x, s2);
%	printf('KB 2d vs nd max%%diff = %g', max_percent_diff(Y2,Y.kb))
	printf('KB 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.kb.p))

	s2 = nufft2_init(om, oldarg{:}, n_shift, 0, 'best');
	printf('tuned 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.mmtu.p))

	s2 = nufft2_init(om, oldarg{:}, n_shift, 0); % minmax
	printf('user 2d vs nd max%%diff = %g', max_percent_diff(s2.p,s.mm.p))

 end
return
end

Nd = [st.xy.Nd st.z.Nd];
Kd = [st.xy.Kd st.z.Kd];

% applying scaling factors
sn = reshape(st.xy.sn(:) * st.z.sn', Nd);
x = x .* sn;
  
% Wz direction interpolation
Xk = reshape(x, [Nd(1)*Nd(2) Nd(3)])';
Xk = fft(Xk, Kd(3), 1); %oversampled FFT, padded at end
if ~isvar('st.z.interp_table')
	Xk = st.z.p * Xk;					% [Mz, Nd(1)* Nd(2)]
else
	Xk = feval(st.z.interp_table, st.z, Xk);
end

% Wx-Wy plane interpolation
Xk = conj(Xk'); %Note "'" is conjugate transpose
Xk = reshape(Xk, [Nd(1:2) st.z.M]);
Xk = fft(Xk, Kd(1), 1); %oversampled FFT, padded at end
Xk = fft(Xk, Kd(2), 2);

nsa = size(st.xy.om,1);
X = zeros(nsa, st.z.M);
for ii = 1 : st.z.M
    if isvar('st.xy.interp_table')
        %load desired frequency locations for different Wx-Wy planes
        st.xy.om = om_xy(nsa*(ii-1)+1:nsa*ii,:); 
        X(:,ii) = feval(st.xy.interp_table, st.xy, col(Xk(:,:,ii)));
    else
        error 'only table-based nufft handled.'
    end
end
clear Xk

