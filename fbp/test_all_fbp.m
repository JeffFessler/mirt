% test_all_fbp

% todo:
% cuboid_im test
% cuboid_proj test

list = {
'cbct_back test'
'ct_geom test'
'image_geom test'
'sino_geom test'
'cylinder_proj test'
'df_example1'
'ellipse_im test'
'ellipse_sino test'
'ellipsoid_proj test'
'ellipsoid_im test'
'fbp_fan_arc_example'
'fbp_fan_arc_point'
'fbp_fan_flat_example'
'fbp_ramp test'
'fbp2_sino_filter test'
'fbp2_example'
'feldkamp_example'
'jaszczak1 test'
'ir_radon_zwart_powell test'
'rebin_helix test' % helix_example
'rect_im test'
'rect_sino test'
%'sphere_proj test'
};

im nan-fail
run_mfile_local(list)
