% block_fatrix_test.m
% Test the block_fatrix object

if 1 || ~isvar('A5'), printm 'setup'
	rng(0)
%	A1 = Gsparse(sparse(rand(10,20)));
	A1 = rand(10,20);
	A2 = magic(20);
	A3 = rand(10,30);
	A4 = Gdft('mask', true(5,4), 'samp', rand(5,4) > 0.5); % for gram
	A5 = rand(size(A4));

	Ac = block_fatrix({A1, A2}, 'type', 'col');
	Ad = block_fatrix({A1, A4}); % diag
	Ak = block_fatrix({A1}, 'type', 'kron', 'Mkron', 2);
	Ar = block_fatrix({A1, A3}, 'type', 'row');
	As = block_fatrix({A4, A5}, 'type', 'sum');
end

if 1, printm 'basic Fatrix tests'
	tester = @(A, mask, name) ...
		Fatrix_test_basic(A, mask, 'complex', 1, 'name', name)
	tester(Ac, true(20,1), 'A:col')
	tester(Ad, true(40,1), 'A:diag')
	tester(Ak, true(40,1), 'A:kron')
	tester(Ar, true(50,1), 'A:row')
	tester(As, true(20,1), 'A:sum')
end

if 1, printm 'adjoint tests'
	tester = @(A) test_adjoint(A, 'complex', 1);
	tester(Ac);
	tester(Ad);
	tester(Ak);
	tester(Ar);
	tester(As);
end

if 1 % test with a "tomo" object
	ig = image_geom('nx', 8, 'ny', 7, 'dx', 1);
	ig.mask = ig.circ > 0;
	sg0 = sino_geom('par', 'nb', 9, 'na', 10, 'dr', 1);
	sg1 = sg0; sg1.na = 6; sg1.orbit = sg0.orbit / sg0.na * sg1.na;
	sg2 = sg1; sg2.na = sg0.na - sg1.na; sg2.orbit_start = sg1.orbit;
		sg2.orbit = sg0.orbit - sg1.orbit;
	tmp = {'square/strip', 'Ltab', [1000], 'strip_width', sg0.dr};
	At0 = Gtomo2_table(sg0, ig, tmp);
	At1 = Gtomo2_table(sg1, ig, tmp);
	At2 = Gtomo2_table(sg2, ig, tmp);
	Atb = block_fatrix({At1, At2}, 'type', 'col', 'tomo', true);

	if 1
		xx = ig.circ;
		yt = cat(2, At1 * xx, At2 * xx);
		yb = Atb * xx;
		jf_equal(yb, yt)
		jf_equal(yb, At0 * xx)
	end

	if 1
		yy = Atb * ig.circ;
		yc = mat2cell(yy, [sg1.nb], [sg1.na sg2.na]);
		xx = At1' * yc{1} + At2' * yc{2};
		x0 = At0' * yy;
		equivs(x0, xx)
		xb = Atb' * yy;
		jf_equal(xb, xx)
	end

%	equivs(full(At0), full(Atb))
	tester_tomo2(Atb, ig.mask, 'A2', At0, 'halt', 0, 'nblock', 0)
end

if 1 % test 'col'
	x = [1:ncol(A1)]';
	y1 = A1 * x;
	y2 = A2 * x;
	yy = Ac * x;
%	printm('col forw error %g', max_percent_diff([y1; y2], yy))
	jf_equal([y1; y2], yy)

	x1 = A1' * y1;
	x2 = A2' * y2;
	xx = Ac' * [y1; y2];
%	printm('col back error %g', max_percent_diff(x1+x2, xx))
	jf_equal(x1+x2, xx)
end

if 1 % test 'diag'
	x1 = [1:ncol(A1)]';
	x2 = [1:ncol(A4)]';
	x = [x1; x2];
	y1 = A1 * x1;
	y2 = A4 * x2;
	yy = Ad * x;
%	printm('diag forw error %g', max_percent_diff([y1; y2], yy))
	jf_equal([y1; y2], yy)

	x1 = A1' * y1;
	x2 = A4' * y2;
	xx = Ad' * yy;
%	printm('diag back error %g', max_percent_diff([x1; x2], xx))
	jf_equal([x1; x2], xx)

	Td = build_gram(Ad, [], 0);
	y1 = Td * xx;
	y2 = [A1' * A1 * x1; A4' * (A4 * x2)];
%	printm('diag gram error %g%%', max_percent_diff(y1, y2))
	jf_equal(y1, y2)
end

if 1 % test 'kron'
	xx = rand(ncol(A1), Ak.arg.Mkron);
	clear xt yt
	for ic=1:Ak.arg.Mkron
		yt(:,ic) = A1 * xx(:,ic);
	end
	yy = Ak * xx(:);
%	printm('kron forw error %g', max_percent_diff(yt(:), yy))
	jf_equal(yt(:), yy)

	yy = reshapee(yy, [], Ak.arg.Mkron);
	for ic=1:Ak.arg.Mkron
		xt(:,ic) = A1' * yy(:,ic);
	end
	xx = Ak' * yy(:);
%	printm('kron back error %g', max_percent_diff(xt(:), xx))
	jf_equal(xt(:), xx)
end

if 1 % test 'row'
	x1 = [1:ncol(A1)]';
	x2 = [1:ncol(A3)]';
	y1 = A1 * x1;
	y2 = A3 * x2;
	yy = Ar * [x1; x2];
%	printm('row forw error %g', max_percent_diff(y1+y2, yy))
	equivs(y1+y2, yy)

	x1 = A1' * yy;
	x2 = A3' * yy;
	xx = Ar' * yy;
%	printm('row back error %g', max_percent_diff([x1; x2], xx))
	equivs([x1; x2], xx)
end

if 1 % test 'sum'
	x = [1:ncol(A4)]';
	y1 = A4 * x;
	y2 = A5 * x;
	yy = As * x;
%	printm('sum forw error %g', max_percent_diff(y1+y2, yy))
	jf_equal(y1+y2, yy)
	x1 = A4' * yy;
	x2 = A5' * yy;
	xx = As' * yy;
%	printm('sum back error %g', max_percent_diff(x1+x2, xx))
	jf_equal(x1+x2, xx)
end
