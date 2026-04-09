%% test tridiag mex
%% check value against \ and ir_apply_tridiag

% todo: compile on ir71 with
% icc -avx2
% and auto vectorize --vec -report=7?
% icc -help
% restrict? --aligned ?
% floating point precision options? -fp-model fast=2
% use simd ?

if 0 % compile if needed, but prefer using Makefile
	mex -O -v CFLAGS="\$CFLAGS -O3 -std=c99 -DMmex -Wall -DNeed_ir_tridiag_inv_mex_gateway" ...
		-I../def/ ir_tridiag_inv_mex.c
	return
end

if 1
	N = 256;
	M = 258;
	scale = 10;
	nthread = int32(jf('ncore'));
	nrep = 4;
	fit = @(x,y) norm(x(:) - y(:)) / norm(x(:)) * 100; % nrmse

	for jj = 1:nrep
		rng(jj);
		% hyperthreading means up to 160 on mpel
		d = scale * randn(N,M);
		d = d + 1i * scale * rand(N,M); % rhs
		a = scale * rand(N-1,1) - scale/2;
		b = scale * rand(N,1) - scale/3;
		c = scale * rand(N-1,1) - scale/2;

		d = single(d);
		a = single(a);
		b = single(b);
		c = single(c);

		T = diag(a,-1) + diag(b) + diag(c,1);
		pr cond(T)
		xb = T \ d;

		if 1
			x1 = ir_apply_tridiag_inv(a, b, c, d);
			pr fit(xb,x1)
		end

		if 1
			try
				x2 = ir_tridiag_inv_mex(a, b, c, d, nthread);
			catch
				printm('ir_tridiag_inv_mex failed');
			end

			pr fit(xb,x2)
		end
	end
end
