% test_all_reg.m

list = {
	'C2sparse test'
	'Cdiff_test'
	'Csparse test'
	'Robject_test'
	'ir_reg_diff_zeroed test'
%	'ir_hct_rgrad'
	'penalty2_design test'
	'penalty2_nuyts test'
%	'penalty_funcs test'
	'potential_func test'
	'potential_shift test'
	'qpwls_psf test'
	'tomo2_beta_test'
	'ugibb_form test'
	'Reg1 test'
	'Cdiff1_test' % slow so last
%	'Cdiff1_tune'
};

run_mfile_local(list)
