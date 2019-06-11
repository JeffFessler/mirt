% edge_response_1d
% show how edge response is nonlinear for edge-preserving regularization

if ~isvar('ytrue'), printm 'setup'
	nx = 2^5;
	ny = nx;
	ny = 20*nx;
	steps = [1 2 4 8];
	nstep = numel(steps);
	xtrue = zeros(nx, ny, 'single');
	ny_per = ny / nstep;
	for ii=1:nstep
		iy{ii} = [1:ny_per] + (ii-1)*ny_per;
		xtrue(nx/2+1:end, iy{ii}) = steps(ii); % for-step
%		xtrue(round(ii/(nstep+1)*nx), iy{ii}) = steps(ii); % impulse
	end

	psf = [1 2 1]';
	mask = true(nx,ny);
	A = Gblur(mask, 'psf', psf);
	ytrue = A * xtrue;

	if 1
		im('notick', xtrue, ' ')
%		ir_savefig edge_response_1d_xtrue
	prompt
	end

	if 1
		im('notick', ytrue, ' ')
%		ir_savefig edge_response_1d_ytrue
	prompt
	end
end

if ~isvar('R'), printm 'R'
	l2b = -0;
	bet = 2^l2b; % essentially product of beta * sigma^2
%	f.pot = {'quad'};
	f.pot = {'hyper3', 1.0};

	% todo: use periodic boundary conditions?
% todo: compare noise with step
	R = Reg1(mask, 'beta', bet, 'type_penal', 'mat', ...
		'offsets', 1, 'pot_arg', f.pot); % 1d regularization
	if 1
		qpwls_psf(A, R, 1, mask);
	prompt
	end
end

rng(0)
%sig = 0.0;
sig = 0.1;
yi = ytrue + sig * randn(size(ytrue));

if ~isvar('xh'), printm 'xh'
	xh = pwls_pcg1(xtrue(mask), A, 1, yi(:), R, 'niter', 30);
	xh = embed(xh, mask);

	if 1
		im('notick', xh, ' ')
%		cbar, colormap(jet)
%		ir_savefig edge_response_1d_xh
	prompt
	end
end

if 1 % plots
	tmp = zeros(nx, nstep);
	stds = zeros(nx, nstep);
	leg = cell(nstep,1);
	for ii=1:nstep
		tmp(:,ii) = mean(xh(:, iy{ii}),2) / steps(ii); % normalize
		stds(:,ii) = std(xh(:, iy{ii}),0,2);
		leg{ii} = num2str(steps(ii));
	end
	if 0 % for-step
		plot(-nx/2:nx/2-1, tmp, '.-')
		axis([0*nx/2+[-1 1]*5 -0.1 1.1])
	else
		plot(tmp, '.-')
		axis([0 nx -0.1 1.1])
		plot(stds, '.-')
	end
	ir_legend(leg)
%	ytick([0 1])
	xlabel 'horizontal location'
	ylabel 'normalized profile'
%	ir_savefig eps_c edge_response_1d_p1
end
