% fig_interp_alf1.m
% visualize the interpolator for optimized alphas
% use for Fig. 7 of Fessler/Sutton 2003 NUFFT paper

if 1 || ~isvar('kerns')
	J = 6;
	M = 2;
	N = ceil(J/M);	% just big enough so K=J
	K = M * N;
	alphas = {[1], [1 -0.46], [1 -0.57 0.14]};
	betas = [0 0.19 0.43];
	gam = 2*pi/K;

	kk = 0;
	ks = linspace(-J/2,J/2, J*100+1)';
	kerns = zeros(length(ks), length(alphas));

	for ii=1:length(alphas)
		alf = alphas{ii};
		beta = betas(ii);

		kern = nufft1_kernel(ks+kk, N, J, K, alf, beta);
		sn = nufft_scale(N, K, alf, beta);

		kerns(:,ii) = kern * mean(sn);
	end
end

%
%	plot compare to dirichlet
%
if 1
	ki = -floor(J/2):floor(J/2);
	clf, plot(...
		ks, kerns(:,1), 'g-', ...
		ks, kerns(:,2), 'c--', ...
		ks, kerns(:,3), 'y-.', ...
		ki, ki==0, 'r.')
%		ks, nufft_diric(ks, N, N), 'm:', ...

	axis tight
	axisy(-0.25, 1.05)
	leg = {};
	text(0.7, 0.98, 'L')
	text(1.5, 0.98, '\alpha')
	for ii=1:length(alphas)
		tmp = sprintf('%g ', alphas{ii});
		tmp = tmp(1:end-1);
		leg{ii} = sprintf('L=%d', ii-1);
		text(0.7, 0.96-0.09*ii, sprintf('%d   (%s)', ii-1, tmp))
	end
	hold on
	plot([0.7 2.7], 0.94*[1 1], 'r-')
	plot(0.95*[1 1], [0.65 1.02], 'r-')
	hold off
%	leg{4} = 'Dirichlet';
	legend(leg, 2)

	xlabel '\omega / \gamma', ylabel 'Interpolator: [R r(\omega)]_1'
	title(sprintf('Equivalent min-max interpolator for J=%d, K/N=%d', J, K/N))

%	printc('fig_interp_alf1')
return
end
