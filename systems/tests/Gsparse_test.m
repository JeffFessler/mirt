% Gsparse_test.m
% test the Gsparse object

if 1
	ig = image_geom('nx', 20, 'ny', 18, 'dx', 1);
	nx = ig.nx;
	ny = ig.ny;
	ig.mask = ig.circ > 0;
	nb = 21; na = 9;
	Al = Glinear('nx', ig.nx, 'ny', ig.ny, 'nb', nb, 'na', na, ...
		'ray_pix', 1, 'mask', ig.mask);
	im plc 3 3; im(1, ig.mask, 'mask')
	if im, im subplot 2, spy(Al), title 'Al', end
end


if 1, printm 'test ordinary sparse matrix'
	A1 = Gsparse(Al, 'mask', ig.mask, 'odim', [nb na]);

	x = single(ig.mask);
	y1 = A1 * x;
	im(4, y1, 'A*x')

	x1 = A1' * y1;
	im(7, x1, 'A''y')

	if isa(A1, 'Fatrix')
		Fatrix_test_basic(A1, ig.mask);
	else
		fatrix2_tests(A1)
	end
	test_adjoint(A1);
	tester_tomo2(A1, ig.mask)
end


if 1, printm 'test multiple rhs'
	j = find(A1.arg.mask);
	j = find(j == sub2ind([nx ny], nx/2, ny/2+2));
	t = reshape(A1(:,j), nb, na);
	im(3, t, 'A(:,j)'), cbar

	tmp = A1(:,[j j+1]);
	tmp = full(tmp); % trick: needed for fatrix2
	t = reshape(tmp, [nb na 2]);
	im(6, t, 'A(:,[j j+1])')

	i = sub2ind([nb na], round(nb/2), round(na/2));
	tmp = A1([i i+2],:)';
	tmp = full(tmp); % trick: needed for fatrix2
	t = embed(tmp, ig.mask);
	im(9, t, 'A([i i+1],:)')
end


if 1, printm 'test blocks'
	nblock = 4;
	Ab = Gblock(A1, nblock);
	x = single(ig.mask);
	istart = 3;
	ia = istart:nblock:na;
	y1 = Ab{istart} * x;
	y2 = Ab * x;
	if any(col(y1-y2(:,ia))), error 'bug', end

	x1 = Ab{istart}' * y1;
	y2 = zeros(size(y2));
	y2(:,ia) = y1;
	x2 = Ab' * y2;
	if any(abs(col(x1-x2)) > 2e-6), error 'bug', end
end


if 1, printm 'test power'
	y1 = A1.^2 * x(ig.mask);
	y2 = Al(:,ig.mask).^2 * double(x(ig.mask));
	if any(abs(col(y1-y2)) > 2e-6), error 'bug', end
end

if 1, printm 'test cell args'
	[i j s] = find(Al);
	s = dsingle(s);
	A2 = Gsparse({i,j,s,nb*na,nx*ny}, 'mask', ig.mask, 'odim', [nb na]);

	x = single(ig.mask);
	y1 = A2 * x;
	y2 = reshape(A2 * x(ig.mask), [nb na]);
	im(5, y1, 'A*x')
	if any(col(y1-y2)), error 'bug', end

	x1 = A2' * y1;
	x2 = embed(A2' * y1(:), ig.mask);
	im(8, x1, 'A''y')
	if any(col(x1-x2)), error 'bug', end

	test_adjoint(A2);
end


if ~has_aspire
	printm 'skipping since no aspire'
	return
end

if 1, printm 'test file.wtf'
	f.dir = test_dir;
	f.wtf = [test_dir 't.wtf'];
	delete(f.wtf)
	wtf_write(f.wtf, Al, nx, ny, nb, na)

	A3 = Gsparse(f.wtf, 'mask', ig.mask);

	x = single(ig.mask);
	y1 = A3 * x;
	y2 = reshape(A3 * x(ig.mask), [nb na]);
	im(6, y1, 'A*x')
	if any(col(y1-y2)), error 'bug', end

	x1 = A3' * y1;
	x2 = embed(A3' * y1(:), ig.mask);
	im(9, x1, 'A''y')
	if any(col(x1-x2)), error 'bug', end

	test_adjoint(A3);
end
