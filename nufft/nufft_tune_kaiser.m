% nufft_tune_kaiser
% tune the shape parameter for kaiser-bessel interpolation in 1D NUFFT
% to minimize the worst-case error

% explore scaling factors
if ~isvar('N') && 0
	N = 2^8; K = 2*N;
	J = 6; alpha = 2.34 * J; kb_m = 0;
	n = [0:(N-1)]'-(N-1)/2;
	kernel = kaiser_bessel('inline', [], alpha, kb_m, K/N);
	sn_zn = reale(1 ./ nufft_interp_zn(0, N, J, K, kernel)); % discrete way
	sn_ft = 1 ./ kaiser_bessel_ft(n/K, J, alpha, kb_m, 1); % cont. FT way
	clf, subplot(121)
	plot(n, sn_zn, 'b.', n, sn_ft, 'r-'), legend('zn', 'FT')
	subplot(122), plot(n, sn_zn ./ sn_ft - 1)
prompt
end

if ~isvar('err.zn')
	N = 2^8;			% large N to reduce effect of N
	n = [0:(N-1)]'-(N-1)/2;
	K_N = 2;			% oversampling factor
	Jlist = [3:16]';
%	Jlist = [2]';			% for m=0, a=2.5 is best
%	mlist = linspace(-2,2,21)';	% J=7 with this list showed m=0
	mlist = 0;
	alist = linspace(2.2,2.4,21)';
	rlist = [0:20]'/40;		% fractions of gamma

	% range of useful alpha values for each J
%	alf_min = sqrt((pi*Jlist/4).^2 - (5.764)^2);
%	alf_min(Jlist < 8) = 0;
%	alf_max = sqrt((pi*Jlist).^2 - (5.764)^2);

	[aa mm] = ndgrid(alist, mlist);
	err.zn = zeros([length(rlist) length(Jlist) numel(aa)]);
	err.ft = zeros([length(rlist) length(Jlist) numel(aa)]);

	for jj=1:length(Jlist)
		J = Jlist(jj);
		printf('J=%d', J)

		K = K_N * N;
		gam = 2*pi/K;
		om = gam * rlist;

		% kaiser-bessel with various shapes
		for ii=1:numel(aa)
			alf = aa(ii) * J;
			kb_m = mm(ii);
			kernel = kaiser_bessel('inline', [], alf, kb_m, K/N);
			kernel_ft = kaiser_bessel_ft('inline', J, alf, kb_m, 1);

			% interpolator worst-case error
		        err.zn(:,jj,ii) = nufft1_error(om, N, J, K, kernel);
		        err.ft(:,jj,ii) = nufft1_error(om, N, J, K, kernel, ...
				kernel_ft);
%				1 ./ kaiser_bessel_ft(n/K, J, alf, mm(ii), 1));
		end

		if 1
			subplot(4,4,jj)
			tmp = reshape(max(err.zn(:,jj,:), [], 1), size(aa));
			semilogy(alist, tmp)
			title(sprintf('J=%d', J)), axis tight
%			prompt
		end
	end
prompt
end

if 1
	emax.zn = max(err.zn, [], 1);			% [M J AO] -> [1 J AO]
	emax.zn = reshape(emax.zn, length(Jlist), numel(aa)); % [J AO]
	[ebest.zn, ibest.zn] = min(emax.zn, [], 2);	% [J AO] -> J
	[ia_best.zn, im_best.zn] = ind2sub(size(aa), ibest.zn);
	if any(ia_best.zn == 1 | ia_best.zn == length(alist))
		warning 'zn end point F'
	end
	if any(im_best.zn == 1 | im_best.zn == length(mlist)) && length(mlist)>1
		warning 'zn end point m'
	end
	abest.zn = alist(ia_best.zn);
	m_best.zn = mlist(im_best.zn);
end

if 1
	emax.ft = max(err.ft, [], 1);		% [M J AO] -> [1 J AO]
	emax.ft = reshape(emax.ft, length(Jlist), numel(aa)); % [J AO]
	[ebest.ft, ibest.ft] = min(emax.ft, [], 2); % [J AO] -> J
	[ia_best.ft, im_best.ft] = ind2sub(size(aa), ibest.ft);
	if any(ia_best.ft == 1 | ibest.ft == length(alist))
		warning 'ft end point F'
	end
	if any(im_best.ft == 1 | im_best.ft == length(mlist)) && length(mlist)>1
		warning 'ft end point m'
	end
	abest.ft = alist(ia_best.ft);
	m_best.ft = mlist(im_best.ft);
end

if 0
	err.m2 = err.zn(:,:,mm == 2);		% just where m=2
	emax.m2 = max(err.m2, [], 1);		% [M J A] -> [1 J A]
	emax.m2 = reshape(emax.m2, length(Jlist), length(alist)); % [J A]
	[ebest.m2, ia_best.m2] = min(emax.m2, [], 2); % [J A] -> J
	if any(ia_best.m2 == 1 | ia_best.m2 == length(alist))
		warning 'm2 end point F'
	end
	abest.m2 = alist(ia_best.m2);
end

% plot vs m (or f).  suprisingly, m=0 seems best!
if 0
	tmp = reshape(emax.zn, [length(Jlist) size(aa)]);	% [J A O]
%	tmp = permute(tmp, [1 3 2]);	% [J O A]
%	tmp = permute(tmp, [3 2 1]);	% [O A J]
	tmp = permute(tmp, [2 3 1]);	% [A O J]
	clf
	semilogy(mlist, 1*tmp(:,:,1), 'y-o'), xlabel m, axis tight
	semilogy(alist, 1*tmp(:,:,1), 'g-^'), xlabel f, axis tight
return
	semilogy(alist, 1*tmp(:,:,1), 'g-^', ...
		alist, 0*tmp(:,:,2), 'y-o')
end

% plot best "fwhm" vs J
clf
if 1
	subplot(131)
	plot(Jlist, abest.zn, 'c-x', Jlist, abest.ft, 'y-o')
	xlabel J, ylabel \alpha, grid, legend('zn', 'ft')
end

% plot ebest vs J
if 1
	subplot(132)
	semilogy(Jlist, ebest.zn, 'c-x', Jlist, ebest.ft, 'y-o')
	xlabel J, ylabel E_{max}, legend('zn', 'ft')
	titlef('Maximum KB error for $K/N=%g$', K/N)

% ir_savefig c 'fig_?'
end

if 1
	subplot(133)
	plot(Jlist, ebest.zn ./ ebest.ft, '-o'), axis tight
	xlabel J, ylabel 'zn / ft'
return
end

% save the tuned parameter to file
% using trick to deal with J=2
if 1
	file = sprintf('kaiser,m=%d', mlist)
	Jlist = [2; Jlist];
	abest.zn = [2.5; abest.zn];
	abest.ft = [2.5; abest.ft];
	save(file, 'Jlist', 'abest')
end
