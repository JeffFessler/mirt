% Cdiff1_test.m
% test Cdiff1 object

% test small size first for adjoint etc.
if 1
	ig = image_geom('nx', 8, 'ny', 6, 'dx', 1);
	list_class = {'fatrix2', 'Fatrix'};
	for ic = 1:numel(list_class)
	 for order = 0:2
		pr order
		arg = {ig.dim, 'order', order, 'class', list_class{ic}};
		if order == 0
			arg = {arg{:}, 'offset', 0};
		else
			arg = {arg{:}, 'offset', [1 2]};
%			arg = {arg{:}, 'offset', [0 1]}; % for 'diff' option
		end

		Cc = Cdiff1(arg{:}, 'type_diff', 'convn');
%		Cd = Cdiff1(arg{:}, 'type_diff', 'diff'); % todo: finish!
		C1 = Cdiff1(arg{:}, 'type_diff', 'for1');
		Cf = Cdiff1(arg{:}, 'type_diff', 'imfilter');
		Ci = Cdiff1(arg{:}, 'type_diff', 'ind');
		Cm = Cdiff1(arg{:}, 'type_diff', 'mex');
		Cp = Cdiff1(arg{:}, 'type_diff', 'circshift');
		Cs = Cdiff1(arg{:}, 'type_diff', 'sparse');
		Cz = Cdiff1(arg{:}, 'type_diff', 'spmat');

		Cc_f = Cc(:,:);
%		Cd_f = Cd(:,:);
		C1_f = C1(:,:);
		Cf_f = Cf(:,:);
		Ci_f = Ci(:,:);
		Cm_f = Cm(:,:);
		Cp_f = Cp(:,:);
		Cs_f = Cs(:,:);
		% Cz already a matrix

		switch list_class{ic}
		case 'Fatrix'
			test_fun = @(C) Fatrix_test_basic(C, true(ig.dim), 'halt', 0);
		case 'fatrix2'
			test_fun = @(C) fatrix2_tests(C);
		otherwise
			fail 'bug'
		end

		test_fun(Cc)
%		test_fun(Cd)
		test_fun(C1)
		test_fun(Cf)
		test_fun(Ci)
		test_fun(Cm)
		test_fun(Cp)
		test_fun(Cs)
%		test_fun(Cz) % Cz is matrix not Fatrix!

		if 1 && order > 0 % trick: Cc,Cf,Cp require Rweights to match
			wt = Rweights(ig.mask, Cc.arg.offset, ...
				'type_wt', 'array', ...
               			'order', order, 'distance_power', 0);

			Ci_f_wt = diag(wt) * Ci_f;

			Cc_f_wt = diag(wt) * Cc_f;
			Cf_f_wt = diag(wt) * Cf_f;
			Cp_f_wt = diag(wt) * Cp_f;

			jf_equal(Ci_f_wt, Cc_f_wt)
			jf_equal(Ci_f_wt, Cf_f_wt)
			jf_equal(Ci_f_wt, Cp_f_wt)

			try
				jf_equal(Ci_f_wt, Cp_f_wt)
			catch
				im plc 2 3
				Cp_f_wt_t = Cp_f_wt';
				Ci_f_wt_t = Ci_f_wt';
				im(1, Ci_f)
				im(4, Cp_f)
				im(2, Ci_f_wt_t)
				im(5, Cp_f_wt_t)
				im(3, Cp_f_wt_t - Ci_f_wt_t)
				im(6, ig.shape(wt))
				fail 'Cp bug'
			end
		else
			jf_equal(Ci_f, Cc_f)
			jf_equal(Ci_f, Cp_f)
		end

%		jf_equal(Ci_f, Cc_f) % see wt'd version above
		jf_equal(Ci_f, C1_f)
%		jf_equal(Ci_f, Cf_f) % see wt'd version above
		jf_equal(Ci_f, Cm_f)
%		jf_equal(Ci_f, Cp_f) % see wt'd version above
		jf_equal(Ci_f, Cs_f)
		jf_equal(Ci_f, Cz)

		% run tests on small cases
		Cdiff1_test1(Cc)
%		Cdiff1_test1(Cd)
		Cdiff1_test1(C1)
		Cdiff1_test1(Cf)
		Cdiff1_test1(Ci)
		Cdiff1_test1(Cm)
		Cdiff1_test1(Cp)
		Cdiff1_test1(Cs)

	 end
	end
end
