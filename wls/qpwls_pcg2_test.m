% qpwls_pcg2_test.m
% test qpwls_pcg2() and examine benefits of circulant preconditioning
% for a simple Toeplitz matrix

if ~isvar('xinit')
	N = 16; nn=1:N;
	xt = ones(N,1); xt(N/2+1) = 2;
	[i1 i2] = ndgrid(nn,nn);
	hx = @(x) (4 + exp(1i*x)) ./ (1 + abs(x)); % kernel of T
	A = hx(i1 - i2);
	T = A' * A;
	clear i1 i2

	% test Fatrix too
	T = Fatrix(size(T), T);

	yy = A * xt;
	bb = A' * yy;

	mask = true(size(xt));
	C = C2sparse('tight', mask, 4);
	C = C(2:N,:); % since 1D
	l2b = 1;
	C = 2^l2b * C;

	H = T.arg + C'*C;
	printm('Cond before %g after %g', cond(T.arg), cond(H))
	xh = H \ bb;
	if im, clf, plot(nn, [xt real(xh) imag(xh)]), end

	xinit = ones(size(xt));
prompt
end

% no preconditioner
niter = 25;
if ~isvar('x1t')
	it = {'niter', niter, 'isave', 'all'};
	M1 = 1;
	x1g = qpwls_pcg1(xinit, A, 1, yy, C, 'precon', M1, it{:});
	x1t = qpwls_pcg2(xinit, T, bb, C, 'precon', M1, it{:});
	equivs(x1g, x1t)
	if im, plot(nn, real(xh), 'yo', nn, real(x1t), '-'), end
prompt
end


% diagonal preconditioner
if ~isvar('x2t')
	M2 = diag(1 ./ diag(H));
	x2g = qpwls_pcg1(xinit, A, 1, yy, C, 'precon', M2, it{:});
	x2t = qpwls_pcg2(xinit, T, bb, C, 'precon', M2, it{:});
	equivs(x2g, x2t)
	if im, plot(nn, real(xh), 'yo', nn, real(x2t), '-'), end
prompt
end

% circulant preconditioner
if ~isvar('M3')
	M3 = qpwls_precon('circ0', {T}, C, mask);
prompt
end

if ~isvar('x3t')
	x3g = qpwls_pcg1(xinit, A, 1, yy, C, 'precon', M3, it{:});
	x3t = qpwls_pcg2(xinit, T, bb, C, 'precon', M3, it{:});
	equivs(x3g, x3t)
	if im, plot(nn, real(xh), 'yo', nn, real(x3t), '-'), end
prompt
end

if ~isvar('x4t') % test stopper
	x4t = qpwls_pcg2(xinit, T, bb, C, 'precon', M3, it{:}, ...
		'stop_threshold', 0.001);
	jf_equal(x3t(:,1:ncol(x4t)), x4t)
end

if 1 && im % plot showing convergence
	ii = 0:niter;
	semilogy(...
		ii, nrms(x1g, xh), 'c.-', ...
		ii, nrms(x1t, xh), 'c-s', ...
		ii, nrms(x2g, xh), 'm.-', ...
		ii, nrms(x2t, xh), 'm-d', ...
		ii, nrms(x3g, xh), 'y.-', ...
		ii, nrms(x3t, xh), 'y-o')
	xlabel 'iteraton', ylabel 'NRMS distance to solution'
	legend('none, g', 'none, t', ...
		'diagonal, g', 'diagonal, t', ...
		'circulant, g', 'circulant, t')
end
