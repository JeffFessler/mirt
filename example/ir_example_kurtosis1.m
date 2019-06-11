% ir_example_kurtosis1
% example of kurtosis for edge-preserving denoising
% Copyright 2012-07-28, Jeff Fessler, University of Michigan

if ~exist('xtrue')
	N = 2^12;
	xtrue = ones(N,1);
	xtrue(end/2+1:end) = 2;
	plot(xtrue)

	sig = 0.2;
	rng(3)
	yi = xtrue + sig * randn(N,1);

	ax = [1 N 0 3];
	plot(yi, '.-'), axis(ax)

	ah = [0 3];
	hh = linspace(0, 3, 121);
	hp = @(y, c) bar(hh, hist(y, hh), c);

	sf = @(xh) text(.1, get(gca, 'ylim') * [0; 0.7], ... 
		sprintf('\\sigma = %3.2f', std(xh(1:round(0.9*N/2)))));

	k1f = @(xh) mean((xh-mean(xh)).^4) / std(xh)^4 - 3;
	kf = @(xh) text(.1, get(gca, 'ylim') * [0; 0.4], ... 
		sprintf('\\gamma_2 = %3.2f', k1f(xh(1:round(0.9*N/2)))));
end

% use local psf to help select beta
if ~isvar('R1'), printm 'R1'
	f.l2b = 2;
	mask = true(N,1);
	R1 = Reg1(mask, 'beta', 2^f.l2b);
%	qpwls_psf(1, R1, 1, mask, 1, 'loop', 0); % choose beta
end

if ~isvar('init'), printm 'init'
	psf = qpwls_psf(1, R1, 1, mask, 1);
	init = conv2(psf, yi, 'same');
end

	nn = N/2 + (-100:102);
%	nn = 1:N;
	ax = [minmax(nn)' 0 3];
	clf, pl = @(i,j) subplot(420+(i-1)*2+j);

	pl(1,1)
	plot(nn, xtrue(nn), 'g.-'), axis(ax), xtick(N/2)
	title 'true signal'

	pl(2,1)
	plot(nn, yi(nn), 'r.-'), axis(ax), xtick(N/2)
	title 'data with additive white gaussian noise'
	pl(2,2)
	hp(yi, 'r'), axisx(ah), title 'Histogram'
	sf(yi), kf(yi)

	pl(4,1)
	plot(nn, init(nn), 'm.-'), axis(ax), xtick(N/2)
	pl(4,2)
	hp(init, 'm'), axisx(ah)
	sf(init), kf(init)

if 0
	pl(1,2)
	ii = N/2+1+[-13:13];
	plot(ii, xtrue(ii), 'g.-', ii, init(ii), 'm.-')
	axis([minmax(ii)' 0.5 2.5])
end


if ~isvar('xpwls1'), printm 'pwls1'
	f.niter1 = 20;
	xpwls1 = pwls_pcg1(init(mask), 1, 1, yi, R1, 'niter', f.niter1);
	plot(xpwls1), axis(ax)
end

	pl(3,1)
	plot(nn, xpwls1(nn), 'y.-'), axis(ax), xtick(N/2)
	title 'PWLS with quadratic regularization'
	pl(3,2)
	hp(xpwls1, 'y'), axisx(ah)
	sf(xpwls1), kf(xpwls1)

if ~isvar('R2'), printm 'R2' % use small delta to show kurtosis
	R2 = Reg1(mask, 'beta', 6*2^f.l2b, 'pot_arg', {'huber', 0.01});
%	qpwls_psf(1, R2, 1, mask, 1);
end


if 0 || ~isvar('xpwls2'), printm 'pwls2'
	f.niter2 = 400;
	xpwls2 = pwls_pcg1(init(mask), 1, 1, yi(:), R2, 'niter', f.niter2);
end

	pl(4,1)
	plot(nn, xpwls2(nn), 'c.-'), axis(ax), xtick(N/2)
	title 'PWLS with TV (strongly edge-preserving) regularization'
	pl(4,2)
	hp(xpwls2, 'c'), axisx(ah)
	axis([0 3 0 300])
	sf(xpwls2), kf(xpwls2)

if 0
	pl(1,2)
	ii = N/2+1+[-13:13];
	plot(ii, xtrue(ii), 'g.-', ii, xpwls1(ii), 'y.-', ii, xpwls2(ii), 'c.-')
	axis([minmax(ii)' 0.5 2.5])
end

%	ir_savefig c ir_example_kurtosis1a
