% interp_table_test.m

hr = @(t, J) (1-abs(t/(J/2))) .* (abs(t) <= J/2);
hi = @(t, J) (1-abs(t/(J/2))).^2 .* (abs(t) <= J/2);

% 1D
if 1, printm '1d'
	L = 10;
	J = 6;
	s = [-J/2*L:J/2*L]'/L;
	hs = hr(s, J) + 1i * hi(s, J);
	if length(hs) ~= J*L+1, error 'size', end

	if im
		clf, subplot(211)
		plot(s, real(hs), 'c.-', s, imag(hs), 'y.-')
	end

	K = 20;
	ck = zeros(K,1);
	ck(2+1) = 1;
	ck = complexify(ck);
	tm = linspace(-2*K, 2*K, 2001)';
	fm = interp1_table_mex(ck, hs, int32(J), int32(L), tm);

	if im
		subplot(212)
		plot(tm, real(fm), 'c.-', tm, imag(fm), 'y.-')
	end
prompt
end

% 2D
if 1, printm '2d'
	L = [2^5 2^4];
	J = [6 4];
	s1 = [-J(1)/2*L(1):J(1)/2*L(1)]'/L(1);
	s2 = [-J(2)/2*L(2):J(2)/2*L(2)]'/L(2);
	h1 = 0*hr(s1, J(1)) + 1i * hi(s1, J(1));
	h2 = 0*hr(s2, J(2)) + 1i * hi(s2, J(2));
%	h1 = complexify(h1);
%	h2 = complexify(h2);

	if im
		clf, subplot(211)
		plot(	s1, real(h1), 'c.-', s1, imag(h1), 'y.-', ...
			s2, real(h2), 'g.-', s2, imag(h2), 'm.-')
	end

	K = [8 12];
	ck = zeros(K);
	ck(0+1, 0+1) = -1i;
%	ck = complexify(ck);
	t1 = linspace(-2*K(1), 2*K(1), 201)';
	t2 = linspace(-2*K(2), 2*K(2), 199)';
	t1 = linspace(0, K(1), 201)';
	t2 = linspace(0, K(2), 199)';
	[tt1 tt2] = ndgrid(t1, t2);
	tm = [tt1(:) tt2(:)];
%	tic
	fm = interp2_table_mex(ck, h1, h2, int32(J), int32(L), tm);
%	toc
	fm = reshape(fm, size(tt1));

	im(121, t1, t2, real(fm), 'real'), cbar
	im(122, t1, t2, imag(fm), 'imag'), cbar
prompt
end


% 3D
if 1, printm '3d'
	L = [2^5 2^4 2^6];
	J = [5 3 4];
	s1 = [-J(1)/2*L(1):J(1)/2*L(1)]'/L(1);
	s2 = [-J(2)/2*L(2):J(2)/2*L(2)]'/L(2);
	s3 = [-J(3)/2*L(3):J(3)/2*L(3)]'/L(3);
	h1 = 1 * hr(s1, J(1)) + 0i * hi(s1, J(1));
	h2 = 1 * hr(s2, J(2)) + 0i * hi(s2, J(2));
	h3 = 1 * hr(s3, J(3)) + 0i * hi(s3, J(3));
	h1 = complexify(h1);
	h2 = complexify(h2);
	h3 = complexify(h3);

	if im
		clf, subplot(211)
		plot(	s1, real(h1), 'c.-', s1, imag(h1), 'y.-', ...
			s2, real(h2), 'g.-', s2, imag(h2), 'm.-', ...
			s3, real(h3), 'y.-', s3, imag(h3), 'w.-')
	end

	K = [8 12 4];
	ck = zeros(K);
	ck(0+1, 0+1) = 1i;
	ck = complexify(ck);
	t1 = linspace(-2*K(1), 2*K(1), 69)';
	t2 = linspace(-2*K(2), 2*K(2), 89)';
	t3 = linspace(-2*K(3), 2*K(3), 9)';
%	t1 = linspace(0, K(1), 201)';
%	t2 = linspace(0, K(2), 199)';
	t3 = [-1 0 1];
	[tt1 tt2 tt3] = ndgrid(t1, t2, t3);
	tm = [tt1(:) tt2(:) tt3(:)];
	tic
	fm = interp3_table_mex(ck, h1, h2, h3, int32(J), int32(L), tm);
	toc
	fm = reshape(fm, size(tt1));

	im(121, t1, t2, real(fm), 'real'), cbar
	im(122, t1, t2, imag(fm), 'imag'), cbar
	if length(t3) == 1 % compare to 2D
		f2 = interp2_table_mex(ck, h1, h2, ...
			int32(J(1:2)), int32(L(1:2)), tm(:,1:2));
		f2 = reshape(f2, size(tt1));
		max_percent_diff(f2, fm)
	end
end
