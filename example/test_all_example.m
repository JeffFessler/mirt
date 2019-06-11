% test_all_example'

list = {
 'cone_beam_ct_example'
 'denoise_threshold_test'
 'ir_ip_unwrap_example'
 'ir_ct_fan_beam_sqs_vs_lalm'
 'l1_regress_example'
 'mri_cs_ist_example'
 'mri_example'
 'mri_example2'
 'mri_example_3d'
 'mri_example_b0'
 'mri_pixel_size_example'
 'mri_sample_density_example'
 'mri_scale_example'
 'mri_sense_demo1'
 'pet_2z_example'
 'pet_transmission_example'
 'radon_example'
 'recon_limited_angle1'
 'restore_ex2'
 'restore_example'
 'restore_l1'
 'restore_sparse'
 'spect_1d_spot_example'
 'spect_3d_example'
% 'unwrap_sps test' % no, this is a routine for unwrap_example
};

run_mfile_local(list)
