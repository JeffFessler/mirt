%	fig_nufft_kspace.m
%	examine accuray of NUFFT (vs DTFT) for k-space MRI trajectory

	N1 = 128; N2 = N1;
	J1 = 6;	J2 = J1;
	K1 = 2*N1; K2 = 2*N2;
	n_shift = [0 0];
%	n_shift = [N1/2 N2/2];

if 0
	kx = [0:(N1-1)]'/N1;
	ky = 0 * kx;
elseif 0
	nturn = 1;
	t = linspace(0,nturn,600)';
	mag = t / nturn;
	phase = 2 * pi * t;
	kx = mag .* cos(phase);
	ky = mag .* sin(phase);
else
	M = 20000;
	t = linspace(0,1,M)';
	w1 = 2*pi*10;
	w2 = 2*pi*3;
	kx = sin(w1 * t) .* cos(w2 * t);
	ky = sin(w1 * t) .* sin(w2 * t);
end
	omega = 2*pi*[kx ky];

if 0
	plot([-1, 1], [0, 0], 'y-', [0 0], [-1, 1], 'y-', ...
		kx, ky, 'c.', kx, ky, 'g-')
	axis square
	text(1.05, 0, 'kx')
	text(0, 1.05, 'ky')
	axis off
	if 0
		xlabel 'kx'
		ylabel 'ky'
		set(gca, 'xtick', [-1 0 1])
		set(gca, 'ytick', [-1 0 1])
		title 'MRI Rosette k-space Trajectory'
	end
%	print('fig_rosette', '-depsc')
end

if 1
	x = [1:N1]' * [1:N2];
	x = shepplogan(N1,N2,1);
	x = x - 0.98 * mean(x(:));

	tic
	Xd = dtft2(x, omega, n_shift, 1);
%	Xd = dtft_mex(omega', x);
	tim.d = toc;

	tic
	Xb = nufft2_bilinear(omega, x, K1, K2, n_shift);
	tim.b = toc;

	tic
	st = nufft2_init(omega, N1, N2, J1, J2, K1, K2, n_shift, 0);
	tim.pre = toc;

	tic
	Xn = nufft2(x, st);
	tim.n = toc

	tic
	s2 = nufft2_init(omega, N1, N2, 2, 2, K1, K2, n_shift, 0);
	tim.pre2 = toc;

	tic
	Xn2 = nufft2(x, s2);
	tim.n2 = toc

	plot(1:M, real(Xd), 'c-', ...
		1:M, real(Xb), 'g-', ...
		1:M, real(Xn), 'y-')
	legend('DTFT', 'Bilinear', 'NUFFT')
	disp(sprintf('NUFFT? max %g%%', max_percent_diff(Xd, Xn)))
	disp(sprintf('NUFFT2 max %g%%', max_percent_diff(Xd, Xn2)))
	disp(sprintf('Bi-lin max %g%%', max_percent_diff(Xd, Xb)))
end
