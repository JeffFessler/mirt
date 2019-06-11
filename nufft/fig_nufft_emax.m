% fig_nufft_emax
% compare max error of gridding interpolators
% for several choices of scaling factors

if 0 || ~isvar('st'), printm 'kernels'
	J = 5;
	% note: optimizing sn with N=2^8, starting with Gaussian, yielded
	% slightly worse error than nufft with KB scaling factors.
	N = 2^8; % for pictures
	N = 2^5; % for optimization
	K = 2 * N;

	gam = 2*pi/K;
	Nmid = (N-1)/2;

	om = linspace(0, gam, 201)';
%	om = gam * [0:199]'/200;
	M = length(om);

	names = {'Uniform', 'Cosine', 'Gaussian', 'KB', 'KB', 'KB:mm', ''};
	types = {'', '', '', '', 'KB:beatty', 'KB:mm', ''};
	for it=1:length(names)
		st{it} = nufft1_build(J, 'N', N, 'K', K, 'om', om, ...
			'sn', names{it}, 'type', types{it});
	end
end


if 1 % show kernel cores
	arg = {};
	kap = st{1}.kap;
	[kap is] = sort(kap(:)); % trick: avoid jumpts in picture
	colors = {'b-', 'g-', 'y-', 'c-', 'r--'};
%	in_list = 1:5;
	in_list = [2 5];
	leg = {};
	for ii=1:length(in_list)
		in = in_list(ii);
		tmp = st{in}.core(is);
		tmp = tmp / max(tmp(:));
		arg = {arg{:}, kap, tmp, colors{in}};
		leg = {leg{:}, [names{in} ' scaling']};
	end
	if 0
	%	arg = {arg{:}, kap, nufft_diric(kap, N, K, 1), 'r-'};
		arg = {arg{:}, kap, nufft_sinc(kap), 'm-'};
		names{end+1} = 'Sinc';
		in_list = [in_list length(st)+1];
	end
	clf
	h = plot(arg{:}, [-J/2 J/2], [0 0], ':', [0 0], [0 1.2], ':');
	h = h(1:end-2);
	axis tight
	axisy([-0.07 1.07])
	h = reshape(h, ncol(kap), length(in_list));
	h = h(1,:);
	legend(h, leg{:})
%	legend(h, names{in_list})
	xlabel 'k'
	title 'LS\_NUFFT interpolation kernel cores'
	text(-2, 0.9, sprintf('J=%d, K/N=%d', J, K/N))
%	ir_savefig(sprintf('fig_nufft_emax_core_%d', J))
%	grid, zoom on
prompt
end


if 0 || ~isvar('err'), printm 'err'
	err = zeros(M, length(st));
	for in=1:length(st)
		err(:,in) = st{in}.err;
	end

	if 0 % examine approximation
		io = 30;
		ft = st{4}.eon(io,:);
		fa = st{4}.approx(io,:);
		nn = st{4}.nn;
		plot(	nn, real(ft), 'b-', nn, imag(ft), 'b:', ...
			nn, real(fa), 'r-', nn, imag(fa), 'r:')
	return
	end
end

if 1, printm 'errors vs om/gam'
	mm = imin(abs(om/gam-0.2));
	x = om(mm)/gam;
	it_list = [1:6];
	it_list = [1:3 5 7];
	ptype = { '-o', '-+', '-^', '-s', '-d', '-x', '-*' };
	arg = {};
	for jj = 1:length(it_list)
		it = it_list(jj);
		arg = {arg{:}, om(mm)/gam, err(mm,it), ptype{it}};
	end
	semilogy(arg{:})
	hold on
	semilogy(om/gam, err(:,it_list), '-')
	hold off
%	axisy(10 .^ [-13 -5])
	axisy(10 .^ [-J -2])
	legend(names{it_list})
	xlabel 'k' %'\omega/\gamma'
	ylabel 'E'
	text(0.05, 10^-2.2, sprintf('J=%d, K/N=%d', J, K/N))
%	ir_savefig(sprintf('fig_nufft_emax_err_%d', J))
	clear mm x it_list ptype arg
prompt
end

if 0, printm 'nrms for random signal';
	rng(0)
	M1 = 200;
	om1 = (rand(M1,1) - 0.5)*2*pi;
	Fm = rand(M1,1) + 1i * rand(M1,1);
	Gd = Gdsft(om1, N);
%	fn = dtft_mex('dtft,adjoint', om1, Fm, N, 1);
	fn = Gd' * Fm;
	Gn_kb = Gnufft({om1, N, J, K, 'kaiser'});
	fn_kb = Gn_kb' * Fm;
	nrms(fn_kb, fn)

	for it=1:length(st)
		ss{it} = nufft1_build(J, 'N', N, 'K', K, 'om', om1, ...
			'sn', names{it}, 'type', types{it});
		fa{it} = ss{it}.slow1(Fm);
		printf('%s\t\&%g %%', names{it}, 100*nrms(fa{it}, fn))
	end
end

% optimize scaling factors
if ~isvar('sn1')
	arg = {J, 'N', N, 'K', K, 'om', om, 'type', 'minmax', 'sn'};
	it = 3; printm(names{it})
	tmp = st{it}.sn(1);
	cost = @(sn) nufft1_emax(arg{:}, [tmp; sn]);
	sn0 = st{it}.sn(2:end); % initial
%	sn1 = fminunc( cost, sn0 );
	[sn1 eopt flag] = fminsearch( cost, sn0 );
%	log10(cost(st{4}.sn))

	st{end} = nufft1_build(arg{:}, [tmp; sn1]);
	err(:,end) = st{end}.err;
	names{end} = 'Optimized';
	log10(max(err))
	fig_nufft_emax
end
