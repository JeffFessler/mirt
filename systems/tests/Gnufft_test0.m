% Gnufft_test0.m
% Basic tests of small Gnufft object

if ~isvar('A') % create Gnufft class object
	omega = linspace(0, 10*2*pi, 101)'; % crude spiral:
	omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);
	if 1 % 2d
		N = [16 14];
		J = [8 6];
		testadj = @(sys) test_adjoint(sys, 'complex', 1);
	else % 3d
		N = [16 14 10];
		J = [8 6 4];
		omega(:,3) = linspace(-pi, pi, size(omega,1));
		testadj = @(sys) test_adjoint(sys, 'big', 1, 'complex', 1);
	end
	K = 2*N;
	mask = true(N);
	mask(:,end) = false;
	cl = 'Fatrix';
	cl = 'fatrix2';
	A = Gnufft(cl, mask, {omega, N, J, K});

	if 1
		fatrix2_tests(A, 'complex', 1, 'tol_gram', 7e-5)
		testadj(A);
	end
end


if ~isvar('W')
	wi = [1:size(omega,1)]';
	W = Gdiag(wi);
end


if 1, printm 'Gnufft gram'
	T = build_gram(A, wi);
	if isa(A, 'Fatrix')
		Fatrix_test_basic(T, mask, 'complex', 1)
	else
		fatrix2_tests(T, 'complex', 1)
	end
	if length(N) == 3
		testadj(T);
	else
		[t0 t1] = test_adjoint(T, 'complex', 1);
		im plc 1 3, im(1, t0), im(2, t1), im(3, t0 - t1')
	end
prompt
end

if length(N) == 2, printm 'T vs A''WA' % slow in 3D
	T2 = T(:,:);
	Af = A(:,:);
	T1 = Af' * diag(wi) * Af;
	max_percent_diff T1 T2
	im plc 1 3, im(1, T1), im(2, T2), im(3, T1 - T2)
%	equivs(y1, y2)
end
