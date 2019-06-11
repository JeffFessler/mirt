% Gdsft_test.m
% Test the Gdsft object

% create Gdsft object
if ~isvar('A'), printm 'setup Gdsft_test'
	N = [32 30];
	M = 201;
	omega = linspace(0, 10*2*pi, M)'; % crude spiral:
	omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);
	pl = 330;
	if im, clf, subplot(pl+1), plot(omega(:,1), omega(:,2), '.'), end
	n_shift = N/2;

	mask = true(N);
	mask(1:5) = false; % stress

	if 1 % non-mex version
		A = Gdsft(omega, N, 'n_shift', n_shift, 'mask', mask, ...
			'use_mex', 0, 'class', 'fatrix2');
		fatrix2_tests(A, 'complex', 1)
		test_adjoint(A, 'complex', 1, 'big', 1);
	end

	classes = {'Fatrix', 'fatrix2'};
	for ic=1:numel(classes)
		cl = classes{ic};
		A = Gdsft(omega, N, 'n_shift', n_shift, 'nthread', 2, ...
			'mask', mask, 'class', cl);

		fatrix2_tests(A, 'complex', 1)
		test_adjoint(A, 'complex', 1, 'big', 1);
	end
end

% test data
if ~isvar('x'), printm 'setup data'
	x = zeros(N);
	x(5:25,10:25) = 1;
	x(15:20,15:20) = 2;
	x(15,5) = 2;
end

% compare forward
if 1
	yt = dtft(x, omega, 'n_shift', n_shift);
	yd = A * [x(mask) x(mask)]; % test with two
	yd = yd(:,1);
	equivs(yt, yd, 'thresh', 1.3e-6)
%	printf('forward max%%diff = %g', max_percent_diff(yt, yd))
end

% compare adjoint
if 1
	xt = mask .* dtft_adj(yt, omega, N, n_shift);
	xd = A' * [yt yt]; % test with two
	xd = xd(:,1);
	xd = embed(xd, mask);
	equivs(xt, xd)
%	printf('back max%%diff = %g', max_percent_diff(xt, xd))
end

if 1, printm 'test adjoint'
	As = Gdsft(omega, [7 8], 'n_shift', n_shift);
	test_adjoint(As, 'complex', 1);
end

if 1 % look at norm for small problem: sqrt(MN) is loose bound
	pr [norm(As) sqrt(prod(size(As))) sqrt(prod(As.idim))]
	Ad = Gdft('mask', true(As.idim));
	pr [norm(Ad) sqrt(prod(size(As))) sqrt(prod(Ad.idim))]
end
