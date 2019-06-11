% test_all_emission.m

list = {
	'eml_curvature test'
	'eml_osem_example'
	'eml_em_test'
	'eml_osem_test'
	'eml_sps_os_test'
	'eml_psca_test'
	'eql_os_emdp_test'
	'eql_sps_os_test'
	'hole_example1'
	'psf_mismatch_example1'
	'psf_mismatch_example2'
};

run_mfile_local(list)
%run_mfile_local(list, 'pause', 1)
