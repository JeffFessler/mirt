% exp_xform_mex_test
% test mex routine exp_xform_mex

rng(0)
N = 500;
M = 6000;
D = 3;
L = 1;
x = randn(N,L) + 1i * randn(N,L);
u = randn(D,N) + 1i * randn(D,N);
v = randn(D,M) + 1i * randn(D,M);
if 1
	x = double(single(x)); % for v7
	u = double(single(u));
	v = double(single(v));
end

cpu etic
y1 = exp(-v.' * u) * x; % v6
cpu etoc 'time matlab:'

cpu etic
y2 = exp_xform_mex(x, u, v);
cpu etoc 'time double mex:'

d = max_percent_diff(y1, y2);
printf('double max %% diff = %g', d)
if d < 1e-12
	printm('exp_xform_mex double appears to be working')
else
	printm('exp_xform_mex double may have a problem?')
end

if ~is_pre_v7
	xs = single(x);
	us = single(u);
	vs = single(v);

	cpu etic
	y3 = exp_xform_mex(xs, us, vs);
	cpu etoc 'time single mex:'

	d = max_percent_diff(y1, y3);
	printf('double max %% diff = %g', d)
	if d < 1e-4 % why so big?
		printm('exp_xform_mex single appears to be working')
	else
		printm('exp_xform_mex single may have a problem?')
	end
end
