function mri_exp_mult_mex
tmp = which('mri_exp_mult_mex');
[mridir, ~] = fileparts(tmp);
if streq(pwd, mridir)
    fail('Do not run Matlab from within the mri/ directory!')
    % because then the path will put mri_exp_mult_mex.m before the mex file!
end
fail('You need to compile mri_exp_mult_mex.%s', mexext)
