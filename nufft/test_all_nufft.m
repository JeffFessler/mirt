 function test_all_nufft
%function test_all_nufft
%|
%| Test all NUFFT routines, including private ones.
%| This must be a function, not a script, to call private functions,
%| although 2016b on a Mac seems to allow call to private from a script.

% test private routines - must be called from this script
try
	nufft_T test
	nufft_interp_zn test
	printm 'private routines passed test'
catch
	fail 'a nufft private routine failed'
end


list = {
	'dtft test'
	'fftn_fast test'
	'ifftn_fast test'
	'interp_table_test'
	'kaiser_bessel_xray test'
%	'my_fftn'
%	'newfft'
	'nufft1_build test'
	'nufft_init test'
	'nufft_scale test'
	'nufft_sinc test'
	'nufft_table_test'
	'nufft test'
};

run_mfile_local(list)
%run_mfile_local(list, 'draw', 1, 'pause', 1)
