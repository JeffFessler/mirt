% dtft_test.m
% test mex versions of DTFT by comparing to m-file versions

for1 = @(om, x) jf_mex('dtft,forward', om', x, int32(1));
adj1 = @(om, Xm, N) jf_mex('dtft,adjoint', om', Xm, int32(N), int32(1));

%
% 1d tests
%
if 1, printm '1d dtft for'
	N1 = 32;
	M = 10;
	om = 2*pi*[0:(M-1)]'/M; % test with uniform frequency locations
	x = [1:N1]';
	Xm = dtft1(x, om);
	Xc = for1(om, x);
	printm('1d forward real mpd = %g%%', max_percent_diff(Xm,Xc))

	X2 = for1([om 0*om], x); % trick: fake as 2d!
	if any(X2 ~= Xc), error bug, end

	x = x - 2i*flipud(x);
	Xm = dtft1(x, om);
	Xc = for1(om, x);
	printm('1d forward cplx mpd = %g%%', max_percent_diff(Xm,Xc))

	X2 = for1(om, [x x]); % L=2 reps
	if any(X2 ~= [Xc Xc]), error bug, end
end

if 1, printm '1d dtft adj'
	xm = dtft_adj(Xm, om, N1, [], 0);
	xc = adj1(om, Xm, N1);
	printm('1d adjoint cplx mpd = %g%%', max_percent_diff(xm,xc))

	x2 = dtft_adj([Xm Xm], om, N1, [], 0);
%	if any(x2 ~= [xm xm]), error bugxm2, end % strangely nonzero!
	x2 = adj1(om, [Xm Xm], N1);
	if any(x2 ~= [xc xc]), error bug, end
end


%
% 2d tests
%
if 1, printm '2d dtft for'
	N1 = 4; N2 = 6;
	x = [[1:N1]'*ones(1,3), ones(N1,N2-3)]; % test signal

	o1 = 2*pi*[0:(N1-1)]'/N1; % test with uniform frequency locations
	o2 = 2*pi*[0:(N2-1)]'/N2; % test with uniform frequency locations
	[o1 o2] = ndgrid(o1, o2);
	omega = [o1(:) o2(:)];
	omega = [omega; [2 7]];

	Xm = dtft2(x, omega, [], 1);
	Xc = for1(omega, x);
	printm('2d forward real mpd = %g%%', max_percent_diff(Xm,Xc))

	x = x - 2i * flipud(x); % complex
	Xm = dtft(x, omega, 'how', 'loop');
	Xc = for1(omega, x);
	printm('2d forward cplx mpd = %g%%', max_percent_diff(Xm,Xc))

	X2 = for1(omega, stackup(x,x));
	if any(X2 ~= [Xc Xc]), error bug, end
end

if 1, printm '2d dtft adj'
	xm = dtft2_adj(Xm, omega, N1, N2, [], 0);
	xc = adj1(omega, Xm, [N1 N2]);
	printm('2d adjoint cplx mpd = %g%%', max_percent_diff(xm,xc))

	if 1 % stack 2d adjoint
		X = [Xm Xm]; % L=2
		x2 = dtft_adj(X, omega, [N1 N2], [], 0);
	%	if any(x2 ~= stackup(xm, xm)), error bug, end % strangely nonzero!
		x2 = adj1(omega, X, [N1 N2]);
		if any(x2 ~= stackup(xc, xc)), error bug, end
	end

	if 0 % real 2d spectrum (not implemented, and probably not needed)
		X = real(Xm);
		xm = dtft2_adj(X, omega, N1, N2, [], 0);
		xc = adj1(omega, X, [N1 N2]);
		printm('2d adjoint real mpd = %g%%', max_percent_diff(xm,xc))
	end
end


%
% 3d tests
%
if 1, printm '3d dtft for'
	N1 = 4; N2 = 6; N3 = 5;
	rng(3)
	x = rand(N1,N2,N3);

	o1 = 2*pi*[0:(N1-1)]'/N1; % test with uniform frequency locations
	o2 = 2*pi*[0:(N2-1)]'/N2;
	o3 = 2*pi*[0:(N3-1)]'/N3;
	[o1 o2 o3] = ndgrid(o1, o2, o3);
	omega = [o1(:) o2(:) o3(:)];
	omega = [omega; [2 7 1]];

	Xm = dtft(x, omega, 'how', 'loop');
	Xc = for1(omega, x);
	printm('3d forward real mpd = %g%%', max_percent_diff(Xm,Xc))

	x = x - 2i * flipdim(x,3); % complex
	Xm = dtft(x, omega, 'how', 'loop');
	Xc = for1(omega, x);
	printm('3d forward cplx mpd = %g%%', max_percent_diff(Xm,Xc))

	X3 = for1(omega, stackup(x,x));
	if any(X3 ~= [Xc Xc]), error bug, end
end

if 1, printm '3d dtft adj'
	xm = dtft_adj(Xm, omega, [N1 N2 N3], [], 0);
	xc = adj1(omega, Xm, [N1 N2 N3]);
	printm('3d adjoint cplx mpd = %g%%', max_percent_diff(xm,xc))

	if 1 % stack 3d adjoint
		X = [Xm Xm]; % L=2
		x2 = dtft_adj(X, omega, [N1 N2 N3], [], 0);
	%	if any(x2 ~= stackup(xm, xm)), error bug, end % strangely nonzero!
		x2 = adj1(omega, X, [N1 N2 N3]);
		if any(x2 ~= stackup(xc, xc)), error bug, end
	end

	if 0 % real 3d spectrum (not implemented, and probably not needed)
		X = real(Xm);
		xm = dtft_adj(X, omega, [N1 N2 N3], [], 0);
		xc = adj1(omega, X, [N1 N2 N3]);
		printm('3d adjoint real mpd = %g%%', max_percent_diff(xm,xc))
	end
end


%
% now test thread speedup
%
if 1 % thread setup
	clear
	N1 = 2^6;
	N2 = 2^5;
	M = 2^12;
	rng(0)
	omega = rand(M, 2);
	s0 = rand(M,1) + 1i * rand(M,1);
	x0 = rand(N1,N2) + 1i * rand(N1,N2);
	nthread = int32(8);

	arg1 = {'dtft,forward', omega', x0};
	arg2 = {'dtft,adjoint', omega', s0, int32([N1 N2]')};
end

if 1, printm 'threads forw'
	tic
	s1 = jf_mex(arg1{:}, int32(1));
	t1 = toc;
	tic
	s2 = jf_mex(arg1{:}, int32(2));
	t2 = toc;

	if any(s1 ~= s2), error bug, end
	printm('thread for time %g %g, speedup = %g', t1, t2, t1 / t2)
	if 1 % vs exact
		se = dtft(x0, omega, 'how', 'loop');
		printm('2d forward exct mpd = %g%%', max_percent_diff(se,s1))
	end
end

if 1, printm 'threads adj'
	tic
	x1 = jf_mex(arg2{:}, int32(1));
	t1 = toc;
	tic
	x2 = jf_mex(arg2{:}, int32(2));
	t2 = toc;
%	printm('2d adjoint mpd2 = %g%%', max_percent_diff(xm,x2))

	if nrms(x2,x1) > 1e-11, warning 'bug-adj-thread'
		max_percent_diff(x2,x1), end
	printm('thread adj time %g %g, speedup = %g', t1, t2, t1 / t2)
	if 1 % vs exact
		xm = dtft2_adj(s0, omega, N1, N2, [], 0);
		printm('2d adjoint exct mpd = %g%%', max_percent_diff(xm,x1))
	end
end
