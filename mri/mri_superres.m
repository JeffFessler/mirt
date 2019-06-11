% mri_superres.m
% NOT DONE!
%
% Explore super-resolution in MRI with rotated propeller sequences
% Inspired by Benoit Desjardins' work.
% Copyright 2003-7-23, Jeff Fessler, The University of Michigan

%
% create k-space samples, in units of 1/mm
%
if ~isvar('kspaces'), printm 'setup kspace'
	fov = 256;	% [mm] typical brain FOV
	N0 = 32;	% nominal image size
	kmax = N0/2*(1/fov);	% for display axes

	k1 = ([-N0/2:N0/2-1]+0.0)/fov;
	[kk1, kk2] = ndgrid(k1, k1);
	kcart = [kk1(:), kk2(:)];

	Nprops = [1 5];	% # of propellers

	for il=1:length(Nprops)
		Nprop = Nprops(il);

		kspace = zeros(Nprop*N0^2, 2);
		for ip=0:Nprop-1
			ang = ip/Nprop * pi;
			rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];
			kspace(ip*N0^2+[1:N0^2],:) = kcart * rot;
			if ip > 0, % avoid repeated DC sampling
				kspace((N0/2+1)^2-1,:) = [];
			end
		end

		kspaces{il} = kspace;
	end, clear ang rot ip il k1 kk1 kk2 kcart

	im plc 2 3
	im('subplot', 1)
	if im
		plot(kspace(:,1), kspace(:,2), '.')
		axis(1.1*[-1 1 -1 1]*kmax), axis square
		xlabel 'k_1 [mm^{-1}]', ylabel 'k_2 [mm^{-1}]'
		title(sprintf('%d k-space samples (%d)', size(kspace,1), Nprop))
	end

	clear kmax kspace
prompt
end


%
% true object
%
if 0 || ~isvar('xtrue'), printm 'setup object'

	% display images with many pixels...
	Ndisp = 256;
	x1d = [-Ndisp/2:Ndisp/2-1] / Ndisp * fov;
	[x1dd x2dd] = ndgrid(x1d, x1d);

	% parameter units all in [mm]
	obj = mri_objects('case1');
	xtrue = obj.image(x1dd, x2dd);

	clim = [0 2];
	im(2, x1d, x1d, xtrue, 'x true', clim), cbar

	clear x1dd x2dd xpar
prompt
end


%
% analytical k-space data
%
if 0 || ~isvar('yd'), printm 'data'

	% noiseless data
	rng(0)
	for ik=1:length(Nprops)
		kspace = kspaces{ik};
		ytrue{ik} = obj.kspace(kspace(:,1), kspace(:,2));
		yd{ik} = ytrue{ik} + 0 * randn(size(ytrue{ik}));
	end
	clear x1dd x2dd xpar kspace ytrue

	[xhatg, yd_g, xg, kg] = mri_grid_linear(kspaces{1}, yd{1}, N0, fov);
	im(3, kg{1}, kg{2}, abs(yd_g), '|y_d|'), cbar

	clear kg yd_g xhatg ik kspace
prompt
end

return

%
% iterative recon
%
if ~isvar('xpcg')

Nlist = [2.^[6:6]];
niter = 20;
beta = 2^-10 * size(kspace,1);
xpcg = {};

for in=1:length(Nlist)
	N = Nlist(in);

	% gridding estimate to initialize iterations
	[xhatg, yhatg, xg] = mri_grid_linear(kspace, yd, N, fov);
	tmp = sprintf('|x| gridding %d', N);
	im(5, xg{1}, xg{2}, abs(xhatg), tmp, clim), cbar

	%
	% create Gnufft class object
	%
	if 1 || ~isvar('G'), printm 'setup Gnufft object'
		omega = 2*pi*kspace*fov/N;
%		minmax(omega(:,1))
		G = Gnufft({omega, [N N], [6 6], 2*[N N], [N/2 N/2]});
	end

	%
	% reconstruct by PCG
	%
	if 1 || ~isvar('xpcg'), printm 'PCG with quadratic penalty'
		mask = true(N);
		R = Robject(mask, 'beta', beta);
		ytmp = yd(:) * (1/fov)^2 * N^2;	% scaling!
		xiter = qpwls_pcg(xhatg(:), G, 1, ytmp, 0, R.C, 1, niter);
		xpcg{in} = embed(xiter(:,end), mask);
		im(6, xg{1}, xg{2}, abs(xpcg{in}), '|x| pcg quad'), cbar
	end
end % fof
end % if

return

if 0
	im plc 2 3
	im(1, x1d, x1d, xtrue, clim, 'x true'), cbar
	xtick([-fov/2 0 fov/2-1])
	ytick([-fov/2 0 fov/2-1])
	for in=1:length(Nlist)
		x1g = [-N/2:N/2-1]'/N * fov;
		tmp = sprintf('N=%d', Nlist(in));
		im(in+1, x1g, x1g, abs(xpcg{in}), tmp, clim), cbar
		xtick([-fov/2 0 fov/2-1])
		ytick([-fov/2 0 fov/2-1])
	end
%	ir_savefig fig_mri_pixel_size_image
return
end


if 1
	clf
	h(1) = plot(x1d, xtrue(:,Ndisp/2+1), 'c-');
	hold on
	llist = {'r:', 'y--', 'g-.', 'r.', 'm-'};
	arg = {'true'};
	for in=1:length(Nlist)
		N = Nlist(in);
		x1g = [-N/2:N/2-1]'/N * fov;
		tmp = xpcg{in};
		h(in+1) = plot(x1g, abs(tmp(:,N/2+1)), llist{in});
		arg = {arg{:}, sprintf('N=%d', N)};
	end
	hold off
	axis([-fov/2 fov/2 0 2.1])
	xlabel 'horizontal position [mm]'
	ylabel '|f(x,0)|'
	set(0, 'DefaultTextFontSize', 12)
%	set(0, 'DefaultAxesFontSize', 14)
	legend(h, arg{:})
%	ir_savefig fig_mri_pixel_size_profile
end
