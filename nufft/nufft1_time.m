%	nufft1_time.m
%	examine time spent in nufft1 sub-routines

if 1
	N1 = 1024 * 128;
	J1 = 6;

	n_shift = 0;
	Nlist = 2 .^ [8];
	Jlist = [8];

	[NN,JJ] = ndgrid(Nlist, Jlist);

	for ii=1:numel(NN)
		N1 = NN(ii);
		J1 = JJ(ii);
		K1 = 2*N1;

		x = [1:N1]';	% test signal
		x = x(:,ones(1,1000));

		gam = 2*pi/K1;
		o1 = linspace(0,gam,N1)';

%profile clear
%profile on -detail operator
%profile on -detail builtin
		tic
		sb = nufft1_init(o1, N1, J1, K1, n_shift, 'best');
		printf('init time = %g', toc)
%profile on -detail builtin
		[Xn, times] = nufft1(x, sb);
%profile report
%profile off

		printf('time fft=%g other=%g overhead=%g', ...
			times(1), times(2), times(2)/times(1)*100)
	end
return
end
