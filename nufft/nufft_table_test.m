% nufft_table_test.m


% 3D test
if 1, printm '3d'
	if 1 || ~isfield(s, 't'), printm 'precompute structure'
%		ktype = 'minmax:kb';
		ktype = 'kaiser';
		Jd = [3 4 3];
		Nd = [10 8 9];
		Ld = [2^8 2^7 2^7];
		Kd = 2 * Nd;
%		n_shift = zeros(size(Nd)); printm 'easy: 0 shift'
		n_shift = [3 7 2]; % stress it

		gam = 2*pi ./ Kd;
		k1 = linspace(-2*Kd(1), 2*Kd(1), 21)';
		k2 = linspace(-2*Kd(2), 2*Kd(2), 18)';
		k3 = linspace(-2*Kd(3), 2*Kd(3), 23)';
		[kk1 kk2 kk3] = ndgrid(k1, k2, k3);
		om = [gam(1)*kk1(:) gam(2)*kk2(:) gam(3)*kk3(:)];

		tic
		s.p = nufft_init(om, Nd, Jd, Kd, n_shift, ktype);
		printm('pre time init %g', toc)

		tic
		s.t = nufft_init(om, Nd, Jd, Kd, n_shift, 'table', Ld, ktype);
		printm('tab time init %g', toc)
	end

	% 3D forward direction
	if 1, printm '3D forward'
		rng(0)
		x = rand(Nd);
		x = zeros(Nd);
		x(4, 2, 3) = 1;

		tic
		Y.d = dtft(x, om, 'n_shift', n_shift); % exact
		printm('dtft3 time %g', toc)

		tic
		Y.p = nufft(x, s.p);
		printm('pre time nufft %g', toc)
		printm('pre max%%diff = %g', max_percent_diff(Y.d, Y.p))

		tic
		Y.t = nufft(x, s.t);
		printm('tab time nufft %g', toc)
		printm('tab max%%diff = %g', max_percent_diff(Y.d, Y.t))
		printm('tab vs pre max%%diff = %g', max_percent_diff(Y.p, Y.t))

		if 1
			im plc 2 3
			im(1, abs(Y.p - Y.d), 'pre - exact'), cbar
			im(2, abs(Y.t - Y.d), 'tab - exact'), cbar
			im(3, abs(Y.t - Y.p), 'tab - pre'), cbar
		end
	end

	% 3D adjoint
	if 1, printm '3D adjoint'
		X = [1:s.t.M]' + 1i * ones(s.t.M,1);

		tic
		c.d = dtft_adj(X, om, Nd, n_shift); % exact
		printm('dtft time adj %g', toc)
		tic
		c.p = nufft_adj(X, s.p);
		printm('pre time nufft adj %g', toc)
		tic
		c.t = nufft_adj(X, s.t);
		printm('tab time nufft adj %g', toc)

		printm('pre max%%diff = %g', max_percent_diff(c.d, c.p))
		printm('tab max%%diff = %g', max_percent_diff(c.d, c.t))
		printm('tab vs pre max%%diff = %g', max_percent_diff(c.p, c.t))

		if 1
			im(4, abs(c.p - c.d), 'pre adj err'), cbar
			im(5, abs(c.t - c.d), 'tab adj err'), cbar
			im(6, abs(c.t - c.p), 'tab-pre adj'), cbar
		end
	end
end


% 2D test
if 1
	if 1 || ~isfield(s, 't'), printm 'precompute structure'
%		ktype = 'minmax:kb';
		ktype = 'kaiser';
		Jd = [4 5];
		Nd = [20 16];
		Ld = [2^12 2^13];
		Kd = 2 * Nd;
%		n_shift = zeros(size(Nd));
		n_shift = [3 7]; % stress it

		gam = 2*pi ./ Kd;
		k1 = linspace(-2*Kd(1), 2*Kd(1), 51)';
		k2 = linspace(-2*Kd(2), 2*Kd(2), 81)';
		[kk1 kk2] = ndgrid(k1, k2);
		om = [gam(1)*kk1(:) gam(2)*kk2(:)];

		tic
		s.p = nufft_init(om, Nd, Jd, Kd, n_shift, ktype);
		printm('pre time init %g', toc)

		tic
		s.t = nufft_init(om, Nd, Jd, Kd, n_shift, 'table', Ld, ktype);
		printm('tab time init %g', toc)
	end

	% forward direction
	if 1, printm 'test forward'
		x = [1:Nd(1)]'*abs([1:Nd(2)]-(Nd(2)/2));

		tic
		Y.d = dtft2(x, om, n_shift); % exact
		printm('dtft2 time %g', toc)

		tic
		Y.p = nufft(x, s.p);
		printm('pre time nufft %g', toc)
		printm('pre max%%diff = %g', max_percent_diff(Y.d, Y.p))

		tic
		Y.t = nufft(x, s.t);
		printm('tab time nufft %g', toc)
		printm('tab max%%diff = %g', max_percent_diff(Y.d, Y.t))
		printm('tab vs pre max%%diff = %g', max_percent_diff(Y.p, Y.t))

		im plc 2 3
		re = @(t) reshape(abs(t), length(k1), length(k2));
		im(1, k1, k2, re(Y.p - Y.d), 'pre'), cbar
		im(2, k1, k2, re(Y.t - Y.d), 'tab'), cbar
		im(3, k1, k2, re(Y.t - Y.p), 'tab-pre'), cbar
	end

	% 2D adjoint
	if 1, printm '2D adjoint'
		X = [1:s.t.M]' + 1i * ones(s.t.M,1);

		tic
		c.d = dtft2_adj(X, om, Nd(1), Nd(2), n_shift); % exact
		printm('dtft2 time adj %g', toc)
		tic
		c.p = nufft_adj(X, s.p);
		printm('pre time nufft adj %g', toc)
		tic
		c.t = nufft_adj(X, s.t);
		printm('tab time nufft adj %g', toc)

		printm('pre max%%diff = %g', max_percent_diff(c.d, c.p))
		printm('tab max%%diff = %g', max_percent_diff(c.d, c.t))
		printm('tab vs pre max%%diff = %g', max_percent_diff(c.p, c.t))

		im(4, abs(c.p - c.d), 'pre adj err'), cbar
		im(5, abs(c.t - c.d), 'tab adj err'), cbar
		im(6, abs(c.t - c.p), 'tab-pre adj'), cbar
	end
end


% 1D test
if 1 || ~isfield(s, 't'), printm '1d'
%	ktype = 'minmax:kb';
	ktype = 'kaiser';
	Jd = 4;
	Nd = 32;
	Ld = 2^14; % table over-sampling
	Kd = 2 * Nd;
	n_shift = zeros(size(Nd));
	% n_shift = Nd/2;

	gam = 2*pi ./ Kd;
	kv = 2*Kd;
	kv = linspace(-kv, kv, 1001)';
	om = gam * kv;

	tic
	s.p = nufft_init(om, Nd, Jd, Kd, n_shift, ktype);
	printm('pre time init %g', toc)

	tic
	s.t = nufft_init(om, Nd, Jd, Kd, n_shift, 'table', Ld, ktype);
	printm('tab time init %g', toc)
end


if 0
	ramp = [1:Nd]';

	Y.dr = dtft(ramp, om, 'n_shift', n_shift);
	Y.de = dtft(eye(Nd), om, 'n_shift', n_shift);
	printm('dtft max%%diff = %g', max_percent_diff(Y.dr,Y.de*ramp))

	Y.pr = nufft(ramp, s.p);
	Y.pe = nufft(eye(Nd), s.p);
	printm('pre max%%diff = %g', max_percent_diff(Y.pr,Y.pe*ramp))

	printm('max%%diff = %g', max_percent_diff(Y.dr,Y.pr))
	printm('max%%diff = %g', max_percent_diff(Y.de,Y.pe))

	e = Y.pr - Y.dr;
	if im
		clf, plot(kv, abs(e), 'c-')
%		clf, plot(kv, sum(abs(e),2), 'c-')
	end
end

% forward direction
if 1, printm '1d forward'
	x = [1:Nd]';
	% x = unitv(Nd, ii);

	Y.d = dtft(x, om, 'n_shift', n_shift); % exact
	cpu etic
	Y.p = nufft(x, s.p);
	cpu etoc 'pre time nufft'
	printm('pre max%%diff = %g', max_percent_diff(Y.d,Y.p))

	tic
	Y.t = nufft(x, s.t);
	printm('tab time nufft %g', toc)
	printm('tab max%%diff = %g', max_percent_diff(Y.d, Y.t))

	if im
		clf
		subplot(211)
		plot(kv, abs(Y.p - Y.d), 'c-'), title 'pre'
		subplot(212)
		plot(kv, abs(Y.t - Y.d), 'y-'), title 'tab'
	end
prompt
end

% adjoint
if 1, printm '1d adj'
	X = [1:s.t.M]' + 1i * ones(size(om));

	c.d = dtft2_adj(X, [om 0*om], Nd(1), 1, [n_shift 0]); % exact
	tic
	c.p = nufft_adj(X, s.p);
	printm('pre time nufft adj %g', toc)
	tic
	c.t = nufft_adj(X, s.t);
	printm('tab time nufft adj %g', toc)

	printm('pre nufft adj vs dtft adj max%%diff = %g', max_percent_diff(c.d,c.p))
	printm('tab nufft adj vs dtft adj max%%diff = %g', max_percent_diff(c.d,c.t))

	if im
		clf
		nn = 0:Nd(1)-1;
		subplot(211)
		plot(nn, abs(c.p - c.d) / mean(abs(c.d)), 'c-'), title 'pre'
		subplot(212)
		plot(nn, abs(c.t - c.d) / mean(abs(c.d)), 'y-'), title 'tab'
	end
end
