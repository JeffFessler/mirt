% Gtomo2_wtmex_test.m
% Test the Gtomo2_wtmex object
% and the Gtomo2_dscmex object
% todo: test both Fatrix and fatrix2 versions

if ~has_aspire
	return
end

if 1, printm 'thread timing of Gtomo2_dscmex with large case'
	igl = image_geom('nx', 512, 'ny', 496, 'dx', 1);
	sgl = sino_geom('par', 'nb', 600, 'na', 100, 'dr', 0.8);
	igl.mask = igl.circ(igl.fov/2) > 0;
	Ald1 = Gtomo2_dscmex(sgl, igl, 'nthread', 1);
	Ald2 = Gtomo2_dscmex(sgl, igl, 'nthread', jf('ncore'));

	if 1 % proj
		x = single(igl.mask);
		cpu etic
		y2 = Ald2 * x;
		t2 = cpu('etoc', sprintf(': %d thread', Ald2.arg.nthread));
		cpu etic
		y1 = Ald1 * x;
		t1 = cpu('etoc', sprintf(': %d thread', Ald1.arg.nthread));
		printm('proj speedup = %g', t1 / t2)
		jf_equal(y1, y2)
	end

	if 1 % back
		y = sgl.ones;
		cpu etic
		x2 = Ald2' * y;
		t2 = cpu('etoc', sprintf(': %d thread', Ald2.arg.nthread));
		cpu etic
		x1 = Ald1' * y;
		t1 = cpu('etoc', sprintf(': %d thread', Ald1.arg.nthread));
		printm('back speedup = %g', t1 / t2)
		jf_equal(x1, x2)
	end
prompt
end

if 1, printm 'thread timing of Gtomo2_wtmex with large case'
	if ~isvar('Alw2')
		if ~isvar('igl')
			igl = image_geom('nx', 512, 'ny', 496, 'dx', 1);
			sgl = sino_geom('par', 'nb', 600, 'na', 100, 'dr', 0.8);
			igl.mask = igl.circ(igl.fov/2) > 0;
		end
		Alw1 = Gtomo2_wtmex(sgl, igl, 'nthread', 1);
		Alw2 = Gtomo2_wtmex(sgl, igl, 'nthread', jf('ncore'));
	end

	if 1 % proj
		x = single(igl.mask);
		cpu etic
		y2 = Alw2 * x;
		t2 = cpu('etoc', sprintf(': %d thread', Alw2.arg.nthread));
		cpu etic
		y1 = Alw1 * x;
		t1 = cpu('etoc', sprintf(': %d thread', Alw1.arg.nthread));
		printm('proj speedup = %g', t1 / t2)
		jf_equal(y1, y2)
	end

	if 1 % back
		y = sgl.ones;
		cpu etic
		x2 = Alw2' * y;
		t2 = cpu('etoc', sprintf(': %d thread', Alw2.arg.nthread));
		cpu etic
		x1 = Alw1' * y;
		t1 = cpu('etoc', sprintf(': %d thread', Alw1.arg.nthread));
		printm('back speedup = %g', t1 / t2)
		equivs(x1, x2, 'thresh', 2e-6) % not equal due to gather
	end
prompt
end

if ~isvar('f.wtf'), printm 'make .wtf'
	ig = image_geom('nx', 22, 'ny', 20, 'dx', 2);
	sg = sino_geom('par', 'nb', 24, 'na', 18, 'dr', 1.8);
	ig.mask = ig.circ(ig.fov/2) > 0;

	f.chat = int32(0);
	f.dsc = [test_dir 't.dsc'];
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');
	f.wtc = strrep(f.dsc, 'dsc', 'wtc');
	arg = aspire_pair(sg, ig, 'support', ig.mask, 'dscfile', f.dsc);
	os_run(sprintf('echo y | wt gen %s row', f.dsc)) % row grouped for OS
	os_run(sprintf('echo y | wt row2col %s %s', f.wtc, f.wtf))
end

if 1, printm 'test wtfmex asp: commands'
	nthread = int32(1);

	bufr = wtfmex('asp:read', f.wtf, f.chat);
	wtfmex('asp:print', bufr)
	tmp = wtfmex('asp:mask', bufr);
	jf_equal(tmp, ig.mask, 'accept_logical_eq_uint8', true)

	mat = wtfmex('asp:load', f.wtf)';
	tmp = wtfmex('asp:load', f.wtc);
	jf_equal(mat, tmp)
	tmp = wtfmex('asp:load', f.wtc, uint8(ig.mask));
	jf_equal(mat(:,ig.mask), tmp)
	tmp = wtfmex('asp:mat', bufr, f.chat)';
	jf_equal(mat, tmp)

	x0 = single(ig.mask);
	y = wtfmex('asp:forw', bufr, nthread, x0, f.chat);
	equivs(y(:), mat * double(x0(:)));

	y0 = single(sg.ones);
	x = wtfmex('asp:back', bufr, nthread, y0, f.chat);
	equivs(x(:), mat' * double(y0(:)));

	x = wtfmex('asp:back2', bufr, nthread, y0, f.chat);
	equivs(x(:), (mat.^2)' * double(y0(:)));

	yb = wtfmex('asp:proj,block', bufr, nthread, x0, int32(1), int32(4), f.chat);
	yt = zeros(size(y), 'single'); yt(:,2:4:end) = y(:,2:4:end);
	jf_equal(yt, yb)

	xt = wtfmex('asp:back', bufr, nthread, yt, f.chat);
	xb = wtfmex('asp:back,block', bufr, nthread, yt, int32(1), int32(4), f.chat);
	jf_equal(xt, xb)

	bufc = wtfmex('asp:read', f.wtc, f.chat);
	x = wtfmex('asp:stayman2', bufc, y0);
	x = wtfmex('asp:nuyts2', bufc, y0);

%	x = wtfmex('asp:pscd', bufc, x?, dqi?, wi?, dj?); % todo pscd

	% now the internally generated version
	tmp = aspire_pair(sg, ig, 'support', 'array');
	bufg = wtfmex('asp:gensys', tmp', 'col', uint8(ig.mask), f.chat);
	wtfmex('asp:print', bufg)
	tmp = wtfmex('asp:mask', bufg);
	jf_equal(tmp, ig.mask, 'accept_logical_eq_uint8', true)

	if 1 % test writing and re-reading
		f.tmp = [test_dir 'tg.wtf'];
		fid = fopen(f.tmp, 'w');
		if (length(bufg) ~= fwrite(fid, bufg)), fail 'fwrite', end
		if (fclose(fid)), fail 'fclose', end
		tmp = wtfmex('asp:read', f.tmp, f.chat);
		jf_equal(bufg, tmp)
	end

	tmp1 = aspire_buff2mat(bufc);
	tmp2 = aspire_buff2mat(bufg);
	jf_equal(tmp1, tmp2)

	tmp = wtfmex('asp:mat', bufg, f.chat);
	jf_equal(mat, tmp)

	yb = wtfmex('asp:forw', bufg, nthread, x0, f.chat);
	equivs(y, yb)
end

if 1, printm 'test Gtomo2_wtmex'
	if ~isvar('Awarg1'), printm 'w: arg mode'
		Awarg1 = Gtomo2_wtmex(arg, 'nthread', 1);
		Awarg2 = Gtomo2_wtmex(arg, 'nthread', jf('ncore'));
		tester_tomo2(Awarg1, ig.mask, 'G2', Awarg2)
		test_adjoint(Awarg1);
		test_adjoint(Awarg2);
	end

	if ~isvar('Awwtf1'), printm 'w: .wtf mode'
		Awwtf1 = Gtomo2_wtmex(f.wtf, 'nthread', 1);
		Awwtf2 = Gtomo2_wtmex(f.wtf, 'nthread', jf('ncore'));
		tester_tomo2(Awwtf1, ig.mask, 'G2', Awwtf2)
		test_adjoint(Awwtf1);
		test_adjoint(Awwtf2);
	end

	if ~isvar('Aw1'), printm 'w: sg,ig mode'
		Aw1 = Gtomo2_wtmex(sg, ig, 'nthread', 1);
		Aw2 = Gtomo2_wtmex(sg, ig, 'nthread', jf('ncore'));
		tester_tomo2(Aw1, ig.mask, 'G2', Aw2)
		test_adjoint(Aw1);
		test_adjoint(Aw2);
	end
prompt
end


if 1, printm 'test Gtomo2_dscmex'
	if ~isvar('Adarg1'), printm 'd: arg mode'
		Adarg1 = Gtomo2_dscmex(arg, 'nthread', 1);
		Adarg2 = Gtomo2_dscmex(arg, 'nthread', jf('ncore'));
		tester_tomo2(Adarg1, ig.mask, 'G2', Adarg2)
		test_adjoint(Adarg1);
		test_adjoint(Adarg2);
	end

	if ~isvar('Adsc1'), printm 'd: .dsc mode'
		Adsc1 = Gtomo2_dscmex(f.dsc, 'nthread', 1);
		Adsc2 = Gtomo2_dscmex(f.dsc, 'nthread', jf('ncore'));
		tester_tomo2(Adsc1, ig.mask, 'G2', Adsc2)
		test_adjoint(Adsc1);
		test_adjoint(Adsc2);
	end

	if ~isvar('Ad1'), printm 'd: sg,ig mode'
		Ad1 = Gtomo2_dscmex(sg, ig, 'nthread', 1);
		Ad2 = Gtomo2_dscmex(sg, ig, 'nthread', jf('ncore'));
		tester_tomo2(Ad1, ig.mask, 'G2', Ad2)
		test_adjoint(Ad1);
		test_adjoint(Ad2);
	end
prompt
end

if 1 && isvar('Aw1') && isvar('Ad1'), printm 'test consistency of dsc with wtf'
	xw = Aw1' * sg.ones;
	xd = Ad1' * sg.ones;
	equivs(xw, xd)

	yw = Aw1 * ig.ones;
	yd = Ad1 * ig.ones;
	equivs(yw, yd)
end
