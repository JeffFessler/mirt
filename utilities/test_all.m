%% test_all.m
%|
%| Attempt to run all test routines
%| to verify completeness of software distribution.
%|
%| I recommend running these tests to verify the completeness
%| of your installation!  If some test fails, please email me
%| the relevant output messages and tell me what version of Matlab
%| and what OS you are using.  Then go ahead and try to use the
%| toolbox because often most things work even if one test fails.
%|
%| If you do not have all the matlab toolboxes (e.g., image, signal, wavelet)
%| then some tests will fail.  But don't worry about it, I am trying to
%| make it work with just plain matlab, so most of it should work even
%| without any extra toolboxes.
%|
%| Jeff Fessler

%double6 double
%disp 'Note: forcing double precision because sparse(double)*single fails :-('

im off-quiet

prompt draw % do all plots without pausing/prompting
%printm 'Note: hit "r" at the prompt to disable subsequent prompts'

diary_file = ['~/' date '-test_all.txt'];
diary(diary_file)

% uncomment the suite(s) you want to test

test_all_mex % see if all the mex files work (this will fail on Windows)

list = {...
'test_all_util', ...
'test_all_graph', ...
'test_all_reg', ...
'test_all_nufft', ...
'test_all_systems', ...
'test_all_emission', ...
'test_all_transmission', ...
'test_all_wls', ...
'test_all_ct', ...
'test_all_mri', ...
'test_all_example'
};

run_mfile_local(list, 'pause', 0)

diary off
