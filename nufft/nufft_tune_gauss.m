% nufft_tune_gauss
% tune the FWHM for gaussian interpolation in 1D NUFFT
% to minimize the worst-case error

if ~isvar('err_gauss')
	N = 2^8;
	K_N = 2;				% oversampling factor
	Jlist = [2:16];
	Plist = linspace(0.2, 0.4, 101);	% sig = parameter * sqrt(J)
	rlist = [0:20]'/40;			% fractions of gamma

	err_gauss.zn = zeros([length(rlist) length(Jlist) length(Plist)]);
	err_gauss.ft = zeros([length(rlist) length(Jlist) length(Plist)]);

	for jj=1:length(Jlist)
		J = Jlist(jj);	printf('J=%d', J)

		K = K_N * N;
		gam = 2*pi/K;
		om = gam * rlist;

		% gaussian with various FWHM
		for ii=1:length(Plist)
			sig = Plist(ii) * sqrt(J);
			[Gfunc, Gfunc_ft] = nufft_gauss('inline', J, sig);

			% interpolator worst-case error
		        err_gauss.ft(:,jj,ii) = nufft1_error(om, N, J, K, ...
				Gfunc, Gfunc_ft);
		        err_gauss.zn(:,jj,ii) = nufft1_error(om, N, J, K, ...
				Gfunc);
		end
	end
end

emax.zn = squeeze( max(err_gauss.zn, [], 1) );		% [M J F] -> [J F]
[emax.zn, ibest.zn] = min(emax.zn, [], 2);		% [J F] -> J
pbest.zn = Plist(ibest.zn);

emax.ft = squeeze( max(err_gauss.ft, [], 1) );		% [M J F] -> [J F]
[emax.ft, ibest.ft] = min(emax.ft, [], 2);		% [J F] -> J
pbest.ft = Plist(ibest.ft);

% plot best sig vs J (on log scale!)
clf
if 1
	log_J = log(Jlist);
	log_sig.zn = log(pbest.zn .* sqrt(Jlist));
	log_sig.ft = log(pbest.ft .* sqrt(Jlist));

	subplot(131)
	plot(log_J, log_sig.zn, 'co', log_J, log_sig.ft, 'yx'), axis tight
	lsline, legend('zn', 'ft')
	xlabel log(J), ylabel 'log(sig best)'
	pol = polyfit(log_J, log_sig.zn, 1)
	printm('\\sigma \\approx %g * J^{%g}', exp(pol(2)), pol(1))
	pol = polyfit(log_J, log_sig.ft, 1)
	printm('\\sigma \\approx %g * J^{%g}', exp(pol(2)), pol(1))
end


% plot emax vs J
if 1
	subplot(132)
	semilogy(Jlist, emax.zn, 'c-o', Jlist, emax.ft, 'y-x')
	xlabel J, ylabel E_{max}
	titlef('Maximum error for $K/N=%g$', K/N)
	legend('zn', 'ft')

	subplot(133)
	plot(emax.zn ./ emax.ft), axis tight, title 'emax: zn / ft'
%	ir_savefig c 'fig_?'
end


% save both sets of tuned sigma's to file
if 1
	Jgauss2 = Jlist;
	Sgauss2.zn = pbest.zn .* sqrt(Jlist); % tuned sig zn
	Sgauss2.ft = pbest.ft .* sqrt(Jlist); % tuned sig ft
	save nufft_gauss2 Jgauss2 Sgauss2
	printm('MOVE nufft_gauss2.mat TO param-data directory!')
end
