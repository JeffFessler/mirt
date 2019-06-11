% fig_error_J.m
% plot maximum normalized error vs J using QR method
% for min-max NUFFT with uniform scaling factors
% This is Fig. 3 in the 2003 IEEE T-SP paper.
% Jeff Fessler, The University of Michigan

if ~isvar('emax')
%	Nlist = 1;
	Nlist = 2^10;
	Mlist = [1.5 2 2.5 3 4 5];	% over-sampling factors
%	Jlist = [2:11];
	Jlist = [2:20];	% neighborhood size

	[JJ, MM, NN] = ndgrid(Jlist, Mlist, Nlist);

	emax = zeros(size(NN));

	alpha = [1]; beta = [];
	for ii=1:numel(NN)
		N = NN(ii);
		J = JJ(ii);
		M = MM(ii);
		printf('M=%d J=%d', M, J)
		K = M * N;
		gam = 2*pi/K;
		om = gam * [0:20]'/40;
		emax(ii) = max(nufft1_err_mm(om, N, J, K, 'qr', alpha, beta));
	end
end


%
% plot min-max errors
%
if 1
	semilogy(Jlist, emax(:,1), '-o', ...
		Jlist, emax(:,2), '-d', ...
		Jlist, emax(:,3), '-^', ...
		Jlist, emax(:,4), '-s', ...
		Jlist, emax(:,5), '-p', ...
		Jlist, emax(:,6), '-+')
%		Jlist(Jlist <= 99), emax(Jlist <= 99, 2), '-d', ...
%		Jlist(Jlist <= 98), emax(Jlist <= 98, 3), '-^', ...
%		Jlist(Jlist <= 96), emax(Jlist <= 96, 4), '-s', ...
%		Jlist(Jlist <= 95), emax(Jlist <= 95, 5), '-p')
%	axis([minmax(Jlist)'+[-0.2 0.2] 1e-5 1e-1])
%	axis tight
%	axisx([minmax(Jlist)'+[-0.2 0.2]])
	axis([minmax(Jlist)'+[-0.2 0.2] 1e-10 1e-0])
	xtick([2:3:20])
	ytick(10.^[-10:2:0])
%	text(2.2, 2e-4, '(for unit-norm signal)')
	xlabel J, ylabel 'E_{max}'
	title(sprintf('Maximum error for \\alpha = (%g)', alpha))
	leg = {};
	for ii=1:length(Mlist)
		leg{ii} = sprintf('K/N=%g', Mlist(ii));
	end
	legend(leg)

%	ir_savefig c 'fig_error_J'
end
