% test_all_ct.m

list = {
'ct_beta_study'
'ct_sys'
%'de_ftab_build'
'de_ftab_curv test'
'de_ftab_fit test'
%'de_ftab_fit_sprad'
'de_ftab_fm test'
%'de_ftab_invert'
%'de_ftab_inv1 test' % 2018-10-01 fails R2018b due to lsqlin solver issue (todo)
%'de_ftab_iwater'
'de_ftab_s_iter test' % slow
%'de_ftab_xform'
'de_ftab test'
%'de_pl_denom'
%'de_pl_obj'
%'de_pl_osps'
'de_poly_eval test'
%'element_density'
%'wls_simplex'
%'xct_poly1_dercurv'
'xray_atten_interp test'
'xray_filters test'
'xray_material_file_name test'
'xray_read_atten test'
'xray_read_dens test'
'xray_read_spectra test'
}

run_mfile_local(list)
