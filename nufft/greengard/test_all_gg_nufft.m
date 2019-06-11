% test_all_gg_nufft.m

list = {
	'ggnuffttest_type1'
	'ggnuffttest_type1_2d'
	'ggnuffttest_type2'
	'ggnuffttest_type2_2d'
	'gg_Gnufft_test'
};

run_mfile_local(list)
%run_mfile_local(list, 'draw', 1, 'pause', 1)
