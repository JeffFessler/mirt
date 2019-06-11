% samsonov_dcf_example.m
%
% Example on reconstruction using Samsonov's method for choosing DCFs.
%
% K. Khalsa, Mar 2006, modified from JF's mri_density_comp.m

fov = [256 256]; % [mm]
nx = 32; ny = nx;
N0 = nx;
if 1
	t = linspace(0, N0/2*2*pi, N0^2+3)';	% crude spiral:
	kspace = N0/2*(1/fov(1))*[cos(t) sin(t)] .* (t(:,[1 1]) / max(t));
	M = size(kspace, 1);
end

if 1
	figure(1)
	im plc 2 3
	im subplot 1
	plot(kspace(:,1), kspace(:,2), '.')
	axis(1.1*[-1 1 -1 1]*N0/2/fov(1)), axis square
	xlabel 'k_1 [mm^{-1}]', ylabel 'k_2 [mm^{-1}]'
	title(sprintf('%d k-space samples', size(kspace,1)))

	% create Gnufft class object
	omega = 2*pi*kspace*fov(1)/N0;
	G = Gnufft({omega, [N0 N0], [6 6], 2*[N0 N0], 1*[N0/2 N0/2], 'kaiser'});
end

% -------------------- %
% Choose DCWs 
% -------------------- %

 
JJ = 4;	% neighborhood
L = 6;

J = Jop({'J', JJ, 'L', L, 'kspace', kspace, 'fov', fov, ...
	'kn.ktype', 'kaiser'});


% Samsonov
wJack = 1 ./ (J * ones(size(J,2), 1));
w0 = wJack;
niter = 100;
[wSamsonov, steps] = samsonov(wJack, niter, J);
wSams = wSamsonov(:,niter);


if 1
	% gridding
%	wi = mri_density_comp(kspace, dtype, G);
	wi = wSams;
	dtype = 'Samsonov';
	im subplot 2
	semilogy(wi, '-r'), title(sprintf('DCF %s', dtype)), hold on
	semilogy(wJack, '-b'), hold off;
	legend('wSams', 'wJack', 'Location', 'Best');

	% true object and analytical k-space data
	obj = mri_objects('case1');

	x1d = [-N0/2:N0/2-1] / N0 * fov(1);
	[x1dd x2dd] = ndgrid(x1d, x1d);
	xtrue = obj.image(x1dd, x2dd);
	ytrue = obj.kspace(kspace(:,1), kspace(:,2));

	clim = [0 2];
	im(4, x1d, x1d, xtrue, 'x true', clim), cbar

	xSams = reshape(1/M * (G' * (wi .* ytrue)), N0, N0);
	xJack = reshape(1/M * (G' * (wJack .* ytrue)), N0, N0);
	xSamsNorm = xSams/ mean2(xSams);
	xJackNorm = xJack / mean2(xJack);
%	im(5, x1d, x1d, real(xSamsNorm), 'xSams, norm'), cbar
%	im(6, x1d, x1d, real(xJackNorm), 'xJack, norm'), cbar

	im(5, x1d, x1d, real(xSams), 'xSams'), cbar
	title(sprintf('xSams, nrms %g%%', 100*nrms(xSams(:), xtrue(:))))
	im(6, x1d, x1d, real(xJack), 'xJack'), cbar
	title(sprintf('xJack, nrms %g%%', 100*nrms(xJack(:), xtrue(:))))
	%sum(xgrid(:)) / sum(xtrue(:))
	printf('nx = %g, M = %g', nx, M)
	printf('xSams nrms %g%%', 100*nrms(xSams(:), xtrue(:)))
	printf('xJack nrms %g%%', 100*nrms(xJack(:), xtrue(:)))
	printf('xSamsNorm nrms %g%%', 100*nrms(xSamsNorm(:), xtrue(:)))
	printf('xJackNorm nrms %g%%', 100*nrms(xJackNorm(:), xtrue(:)))
	printf('nrms btn wJack and wSams %g%%', 100*nrms(wJack, wSams));
	
	im subplot 3
	plot(x1d, xtrue(:,end/2), 'c-', x1d, real(xSams(:,end/2)), '-r', ...
		x1d, real(xJack(:, end/2)), 'b-')
	axis tight
	legend('true', 'Sams', 'Jack', 'Location', 'Best')
	
	
	% nrms error vs iteration %
	sams_nrms = zeros(niter,1);
	for i=1:niter, wi = wSamsonov(:,i);
		xSams = reshape(1/M * (G' * (wi .* ytrue)), nx, ny);
		sams_nrms(i) = nrms(xSams(:), xtrue(:));
	end
	figure(2)
	plot(1:niter, sams_nrms), title('nrms of xSams vs iteration')
	xlabel('iteration number'), ylabel('nrms')
	
end
return
