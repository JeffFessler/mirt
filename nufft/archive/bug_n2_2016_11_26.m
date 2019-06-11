% bug_n2_2016_11_26.m
% debug issue when N=2 and J>K

if 0
	N = [5]; % works
	N = [4]; % rank deficient warning: nufft_alpha_kb_fit (line 47)
	N = [3]; % works
	N = [2]; % works now because nufft_init forces J to be min(K, J); prec. warn
	Gnufft(true([N 1]), {[0], N, [6], 2*N, N/2, ...
		'table', 2^12, 'minmax:kb'});
end

if 0 % 2D
	N = [1 1] * 3; % works
	N = [1 1] * 2;
	Gnufft(true(N), {[0 0], N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'});
end

if 1 % 3D
	N = [20 20 2];
	Gnufft(true(N), {[0 0 0], N, [6 6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'});
end
