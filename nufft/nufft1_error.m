 function [err, sn, kernel] = nufft1_error(om, N1, J1, K1, kernel, sn)
%function [err, sn, kernel] = nufft1_error(om, N1, J1, K1, kernel, sn)
%|
%| Compute worst-case error for each input frequency for 1D NUFFT
%| using specified (inline) `ad hoc' kernel function (e.g., gaussian).
%| This is worst-case for a unit-norm signal of length N1.
%| in
%|	om	[M,1]	digital frequency omega in radians
%|	N1		signal length
%|	J1		# of neighbors used per frequency location
%|	K1		FFT size (should be > N1)
%|	kernel		inline kernel function, args (k,J)
%|				(or choose from built-in examples - see code!)
%|	sn		optional scaling factors (otherwise do-no-harm)
%| out
%|	err	[M,1]	worst-case error over unit-norm signals
%|	sn	[N,1]	scaling factors
%|
%| examples for kernel:
%| linear:		inline('(1 - abs(k/(J/2))) .* (abs(k) < J/2)', 'k', 'J')
%| truncated diric:	inline('sinc(k) .* (abs(k) < J/2)', 'k', 'J')
%|
%| Copyright 2001-12-7, Jeff Fessler, The University of Michigan

% if no arguments, give an example, comparing triangular kernel to min-max
if nargin < 4
	help(mfilename)
	N = 2^7; K = 2*N; gam = 2*pi/K;
	Jlist = [2:10]';
	om = gam * linspace(0,1,101);

	err.linear	= zeros(size(Jlist));
	err.minmaxu	= zeros(size(Jlist));
	err.minmax2	= zeros(size(Jlist));
	err.minmaxo	= zeros(size(Jlist));
	err.minmaxk	= zeros(size(Jlist));	% kaiser sn's
	err.gauss_zn	= zeros(size(Jlist));
	err.gauss_ft	= zeros(size(Jlist));
	err.kaiser	= zeros(size(Jlist));

	for ii=1:length(Jlist)
		J = Jlist(ii);
		printf('J=%d', J)
		err.linear(ii) = ...
			max(nufft1_error(om, N, J, K, 'linear'));
		err.minmaxu(ii) = ...
			max(nufft1_error(om, N, J, K, 'minmax,uniform'));
		err.minmax2(ii) = ...
			max(nufft1_error(om, N, J, K, 'minmax,best,L=2'));
		err.minmaxo(ii) = ...
			max(nufft1_error(om, N, J, K, 'minmax,best'));
		err.gauss_zn(ii) = ...
			max(nufft1_error(om, N, J, K, 'gauss'));
		err.gauss_ft(ii) = ...
			max(nufft1_error(om, N, J, K, 'gauss', 'ft'));
		[tmp, sn] = nufft1_error(om, N, J, K, 'kaiser', 'ft');
		err.kaiser(ii) = max(tmp);
		err.minmaxk(ii) = ...
			max(nufft1_err_mm(om, N, J, K, 'qr', sn));
	end
	clf, semilogy(Jlist, err.linear, 'g-x', ...
		Jlist, err.minmaxu, 'c-+', ...
		Jlist, err.gauss_zn, 'b-*', ...
		Jlist, err.gauss_ft, 'b-o', ...
		Jlist, err.minmax2, 'r-^', ...
		Jlist, err.kaiser, 'm->', ...
		Jlist, err.minmaxo, 'y-<', ...
		Jlist, err.minmaxk, 'w-o'), axis tight
	xlabel J, ylabel 'worst-case error'
	legend('linear', 'min-max, uniform', ...
		'gaussian (zn)', 'gaussian (FT)', ...
		'min-max, best L=2', 'kaiser', 'min-max, optimized', ...
		'min-max, kaiser s', 3)
	clear err
return
end


%
% kernel selection
%
if ischar(kernel)

	% min-max interpolators (uniform, best, etc., see nufft1_err_mm.m)
	if strncmp(kernel, 'minmax,', 7)
		type = kernel(8:end);	% uniform or best or best,L=2 etc.
		[err, sn] = nufft1_err_mm(om, N1, J1, K1, 'qr', type);
		return

	% cos^3-tapered dirichlet
	elseif streq(kernel, 'cos3diric')
		kernel = 'diric(2*pi*k/J, J) .* cos((2*pi*k/J)/2).^3';

	% Dirichlet (truncated)
	elseif streq(kernel, 'diric')
		kernel = 'nufft_diric(k,%d,%d,1) .* (abs(k) < J/2)';
		kernel = sprintf(kernel, N1, N1);

	% gaussian (truncated) with previously numerically-optimized width
	elseif streq(kernel, 'gauss')
		if isvar('sn') && ischar(sn) && streq(sn, 'ft')
			stype = 'ft';
		else
			stype = 'zn';
		end
		[dummy, kernel, kernel_ft] = nufft_best_gauss(J1, K1/N1, stype);
		kernel_ft = inline(kernel_ft, 't');

	% kaiser-bessel with previously numerically-optimized shape
	elseif streq(kernel, 'kaiser')
		[kernel, kb_a, kb_m] = kaiser_bessel('string', J1, 'best', 0, K1/N1);
		kernel_ft = kaiser_bessel_ft('inline', J1, kb_a, kb_m, 1);

	% linear interpolation via triangular function (possibly "wide"!)
	elseif streq(kernel, 'linear')
		kernel = '(1 - abs(k/(J/2))) .* (abs(k) < J/2)';

	else
		fail('unknown kernel "%s"', kernel)
	end

	kernel = inline(kernel, 'k', 'J');

elseif ~streq('inline', class(kernel))
	error 'need inline kernel'
end

gam = 2*pi/K1;

if 0
	%	plot interpolator
	k = linspace(-J/2-1,J/2+1,101);
	clf, subplot(221), plot(k, kernel(k, J))
	xlabel k, ylabel kernel(k), axis tight, grid
end

%
% Compute scaling factors using the "do no harm" strategy.
% This may not be optimal; analytical FT could also be reasonable.
%
if ~isvar('sn')
	sn = 1 ./ nufft_interp_zn(0, N1, J1, K1, kernel);	% [N]
%	sn = 1 ./ mean(nufft_interp_zn([0 1/2], N1, J1, K1, kernel), 2); % alt

elseif isa(sn, 'inline')
	n = [0:(N1-1)]'-(N1-1)/2;
	sn = 1 ./ sn(n/K1);		% [N]

% trick to use Gaussian FT scaling factors
elseif ischar(sn) && streq(sn, 'ft') && isvar('kernel_ft')
	Nmid = (N1-1)/2; % todo!
	n = [0:(N1-1)]' - Nmid;
	sn = 1 ./ kernel_ft(n/K1);		% [N]
	if 0 && J1 > 2
		sn_zn = 1 ./ nufft_interp_zn(0, N1, J1, K1, kernel);	% [N]
		clf, plot(n, [sn reale(sn_zn)])
		keyboard
	end

% trick to use Kaiser-Bessel FT scaling factors
elseif isstruct(sn) && streq(sn.type, 'kaiser')
	Nmid = (N1-1)/2; % todo!
	n = [0:(N1-1)]' - Nmid;
	sn = 1 ./ kaiser_bessel_ft(n/K1, J1, sn.alpha, sn.m, 1);

else
	error 'unsupport scaling factors type'
end

%
% interpolator worst-case error for each frequency (scaled by 1/sqrt(N)),
% from equations (46)-(47) in Fessler&Sutton NUFFT paper, T-SP, Feb. 2003
%
zn = nufft_interp_zn(om/gam, N1, J1, K1, kernel);		% [N,M]
err = sqrt(mean(abs(spdiag(sn, 'nowarn') * zn - 1).^2, 1))';	% [M]
