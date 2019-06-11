% mri_exp_mult_mex_test
% test mri_exp_mult_mex

L = 10;
N = 30;
M = 20;
rng(0)
A = randn(N,L) + 1i * randn(N,L);
ur = randn(N,1);
ui = randn(N,1);
vr = randn(M,1);
vi = randn(M,1);
u = ur + 1i * ui;
v = vr + 1i * vi;

d1 = A' * exp(-ur * v.');
d2 = mri_exp_mult_mex(A, ur, v);
printf('double u real case error: %g%%', max_percent_diff(d1, d2))

d1 = A' * exp(-u * vr.');
d2 = mri_exp_mult_mex(A, u, vr);
printf('double v real case error: %g%%', max_percent_diff(d1, d2))

if 0 % test both complex, both real
	mri_exp_mult_mex(A, ur, vr);
	mri_exp_mult_mex(A, u, v);
end

if ~is_pre_v7
	u = single(u);
	v = single(v);
	ur = single(ur);
	vr = single(vr);
	A = single(A);

	d1 = A' * exp(-ur * v.');
	d2 = mri_exp_mult_mex(A, ur, v);
	printf('single u real case error: %g%%', max_percent_diff(d1, d2))

	d1 = A' * exp(-u * vr.');
	d2 = mri_exp_mult_mex(A, u, vr);
	printf('single v real case error: %g%%', max_percent_diff(d1, d2))
end
