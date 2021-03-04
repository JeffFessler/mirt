% ir_get_data_test.m
%
% This script tests the ir_get_data function by deleting the data in the
% downloads folder and calling the scripts that will re-download the data.
%
% Melissa Haskell, mhask@umich.edu, 2021
%
%

%% First delete data that is in downloads folder

warning('WARNING!! Running this test script will PERMANENTLY delete all data in the data/dowloads folder!')
disp('Only proceed if you want to do this!')
disp('Press any key to continue or CTRL-C to quit.')
pause

% make sure script is being run from correct directory
if ~strcmp(pwd,fileparts(which('ir_get_data_test')))
    error('Error: Must run ir_get_data_test.m from data directory to avoid deleting wrong files.')
else
    try
        rmdir downloads s
    catch
        disp('No files to delete.')
    end
end


%% run scripts that need downloaded data

% example directory scripts
ct_fan_beam_example; clear
pet_transmission_example; clear
mri_example_epi_b0; clear

% mri directory scripts
mri_field_map_reg('test'); clear
mri_phase_denoise('test'); clear
fmap_est_pcg_ls('test'); clear
fmap_est_qm('test'); clear

fprintf('\n\n ir_get_data() tests completed.\n')


%% scripts that require dowloaded data, but currently won't run bc of mex
%  files, so left out of March 2021 version of tests

%%% from example dir
% ir_ct_fan_beam_sqs_vs_lalm; clear
% ir_ct_roi_split1; clear

%%% from mri dir
% mri_field_map_reg_2d('test'); clear


