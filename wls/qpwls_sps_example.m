% qpwls_sps_example.m
%
% Complete example m-file illustrating QPWLS-SPS algorithm and its relatives.
%
% Copyright 1999-4, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), disp 'setup qpwls_sps_example'
	em_wls_test_setup
%	wi = ones(size(wi)); warning 'uniform wi' % test circulant precon
	W = diag_sp(wi(:));
prompt
end


if 0 && ~isvar('vb'), disp 'make vb'
	vb = G' * (W * yi(:));
	im(7, ig.embed(vb), 'vb'), cbar
prompt
end


%
% regularization matrix
%
if ~isvar('R'), disp 'make R'
	kappa = sqrt( (G' * wi(:)) ./ (G' * ones(size(wi(:)))) );
	kappa = ig.embed(kappa);
	im(8, kappa, 'kappa'), cbar

	f.l2b = 9;
	R = Robject(kappa, 'type_denom', 'matlab', 'beta', 2^f.l2b);
prompt
end


%
% SPS iterations
%
if ~isvar('xsps'), disp 'do sps'
	f.niter = 20;
%	xinit = ones(size(xtrue));	% uniform initial image
	xinit = max(xfbp,0);

	xsps = pwls_sps_os(xinit(ig.mask), yi, wi, G, R, ...
		f.niter, inf, [], [], 1);
	xsps = ig.embed(xsps);
	im clf, im(xsps, 'xsps')
prompt
end


%
% OS-SPS iterations
%
if ~isvar('xos'), disp 'do os'
	nsubset_list = [2 4 8 16];
	for is = 1:length(nsubset_list)
		Gb = Gblock(G, nsubset_list(is));
		tmp = pwls_sps_os(xinit(ig.mask), yi, wi, Gb, R, ...
			f.niter, inf, [], [], 1);
		xos{is} = ig.embed(tmp);
		im clf, im(xos{is})
	end
prompt
end


%
% Ordinary CG iterations
%
if ~isvar('xcg'), disp 'xcg'
	xcg = qpwls_pcg(xinit(ig.mask), G, W, yi, 0, ...
			R.C, 1, f.niter, ig.mask);
	xcg = ig.embed(xcg);
	im clf, im(xcg, 'CG')
prompt
end


%
% Preconditioned CG iterations
%
if ~isvar('xpcg'), disp 'xpcg'
	xpcg = qpwls_pcg(xinit(ig.mask), G, W, yi, 0, ...
			R.C, 'circ0', f.niter, ig.mask);
	xpcg = ig.embed(xpcg);
	im clf, im(xpcg, 'PCG')
prompt
end


%
% compare CG and PCG
%
if 1 && im
	cost.sps	= pwls_cost(xsps,	G, W, yi(:), R, ig.mask);
	cost.cg		= pwls_cost(xcg,	G, W, yi(:), R, ig.mask);
	cost.pcg	= pwls_cost(xpcg,	G, W, yi(:), R, ig.mask);
	ii = 0:f.niter-1;
	im clf
	plot(ii, cost.sps, 'y-o', ii, cost.cg, 'g-x', ii, cost.pcg, 'c-+')
	title 'caution: sps uses nonnegativity but cg does not...'
	legend('SPS', 'CG', 'PCG (circ0)')
	axis([0 10 0 2e5])
%prompt
end

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% steepest descent
%
if ~isvar('xsd')
	xsd = ig.zeros(f.niter);
	x = xinit(ig.mask);
	xsd(:,:,1) = ig.embed(x);

	for ii=2:f.niter
		if (rem(ii,10) == 0), disp(sprintf('SD Iteration %d', ii)), end

		grad = G' * (W * (G*x)) + C' * (C * x) - vb; % grad = Hx-b
		ddir = -grad;		% gradient direction
		t = G * ddir;		% fix: can be made more efficient!
		stepsize = -(ddir' * grad) / (t' * W * t + norm(C*ddir).^2);
%		disp(sprintf('ii=%02d step=%g', ii, stepsize))
		x = x + stepsize * ddir;

		xsd(:,:,ii) = ig.embed(x);
	end
	im clf, im(xsd)
return
end

%
% experiment with OS-based preconditioner
%
if ~isvar('xpre')
	nsubset = 8;
	starts = subset_start(nsubset);
	xpre = ig.zeros(f.niter);

	x = xinit(ig.mask);
	xpre(:,:,1) = ig.embed(x);
	denoms = pwls_sps_os([]);
	im(ig.embed(denoms), 'Denoms')

	for ii=2:f.niter
		if (rem(ii,10) == 0), disp(sprintf('Pre Iteration %d', ii)), end

		% even iterations, do SPS-OS
		if (rem(ii,2) == 0)
			x = pwls_sps_os([]);
		% odd iterations, do PSD
		else
			% grad = Hx-b
			grad = G' * (W * (G*x)) + C' * (C * x) - vb;
			if 1
				ddir = xpre(:,:,ii-2);
				ddir = x - ddir(ig.mask);	% OS-based search direction
			else
				ddir = -grad;		% gradient direction
			end

			ddir = 0.9 * ddir/norm(ddir) + 0.1 * (-grad) / norm(grad);

			% stepsize
			t = G * ddir;
			stepsize = -(ddir' * grad) / (t' * W * t + norm(C*ddir).^2);
%			disp(sprintf('ii=%02d step=%g', ii, stepsize))
			x = x + stepsize * ddir;
		end

		xpre(:,:,ii) = ig.embed(x);
	end
	im clf, im(xpre)
return
end


%
% new idea
%
if ~isvar('xnew')
	nsubset = 4;
	starts = subset_start(nsubset);
	xnew = ig.zeros(f.niter);
	x = xinit(ig.mask);
x = xpcg(:,:,end);	x = x(ig.mask);	disp('cheat')
	xnew(:,:,1) = ig.embed(x);
	denoms = pwls_sps_os([]);
	im(ig.embed(denoms), 'Denoms')

	%
	% compute 'partial gradients'
	%
	rset = zeros(length(x), nsubset);		% space for H_k x - b_k
	proj = G * x;
	ctcx = C' * (C*x);

	for ik=1:nsubset
		ia = starts(ik):nsubset:sg.na;
		ii = col(outer_sum(1:sg.nb,(ia-1)*sg.nb));
		Gti = G(ii,:)';
		rset(:,ik) = Gti * col(wi(:,ia) .* (yy(:,ia) - proj(:,ia)));
		rset(:,ik) = nsubset * rset(:,ik) + ctcx;
	end

	denom = sum(denoms,2);
	for ii=2:f.niter
		if (rem(ii,10) == 0), disp(sprintf('new Iteration %d', ii)), end

disp([range(sum(rset,2),2)' norm(sum(rset,2))])
		for ik=1:nsubset
			ia = starts(ik):nsubset:sg.na;
			iis = col(outer_sum(1:sg.nb,(ia-1)*sg.nb));
			Gti = G(iis,:)';

			if 0	% SPS-OS
				t = col(yy(:,ia)) - (x' * Gti)';
				t = Gti * (col(wi(:,ia)) .* t);
				t = nsubset * t - (C' * (C*x));
				delta = t ./ denoms(:,ik);
			end

%			delta = sum(rset,2) ./ denoms(:,ik);
			delta = sum(rset,2) ./ denom;
%			delta = rset(:,ik) ./ denoms(:,ik);
			x = x + delta;

			if 1
				t = Gti * (col(wi(:,ia)) .* (delta' * Gti)');
				t = nsubset * t + (C' * (C*delta));
				rset(:,ik) = rset(:,ik) - t;
			end
		end
		xnew(:,:,ii) = ig.embed(x);
	end
	im clf, im(xnew)
return
end

% plot cost function
if 1
	cost.sps = pwls_cost(xsps,	G, W, yy(:), R, ig.mask);
	cost.sd	= pwls_cost(xsd,	G, W, yy(:), R, ig.mask);
	cost.pcg = pwls_cost(xpcg,	G, W, yy(:), R, ig.mask);
	cost.os2 = pwls_cost(xos.s2,	G, W, yy(:), R, ig.mask);
	cost.os4 = pwls_cost(xos.s4,	G, W, yy(:), R, ig.mask);
	cost.os8 = pwls_cost(xos.s8,	G, W, yy(:), R, ig.mask);
	cost.os16 = pwls_cost(xos.s16,	G, W, yy(:), R, ig.mask);
	cost.pre = pwls_cost(xpre,	G, W, yy(:), R, ig.mask);
	cost.new = pwls_cost(xnew,	G, W, yy(:), R, ig.mask);
	ii = 0:f.niter-1;

	% transform so that SPS is a straight line
	if 1
		cost.pcg	= interp1x(cost.sps, f.niter:-1:1, cost.pcg, 1);
		cost.sd		= interp1x(cost.sps, f.niter:-1:1, cost.sd, 1);
		cost.os2	= interp1x(cost.sps, f.niter:-1:1, cost.os2, 1);
		cost.os4	= interp1x(cost.sps, f.niter:-1:1, cost.os4, 1);
		cost.os8	= interp1x(cost.sps, f.niter:-1:1, cost.os8, 1);
		cost.os16	= interp1x(cost.sps, f.niter:-1:1, cost.os16, 1);
		cost.pre	= interp1x(cost.sps, f.niter:-1:1, cost.pre, 1);
		cost.new	= interp1x(cost.sps, f.niter:-1:1, cost.new, 1);
		% this one must be last!
		cost.sps	= interp1x(cost.sps, f.niter:-1:1, cost.sps, 1);
	end

	close all, plot(ii, cost.sps, 'r-*', ...
		ii, cost.new, 'y--*', ...
		ii, cost.sd, 'g--x', ...
		ii, cost.pcg, 'r--x', ...
		ii, cost.os2, 'g-^', ...
		ii, cost.os4, 'y->', ...
		ii, cost.os8, 'w-<', ...
		ii, cost.os16, 'c-v', ...
		ii, cost.pre, 'm--+')
	xlabel('iteration'), ylabel('transformed cost function')
%	axisy(800, 1200)
	axisy(-10, 75)
	axisx(0, 30)
	legend('SPS', 'NEW', 'SD', 'PCG-I', 'OS2', 'OS4', 'OS8', 'OS16', 'PRE', 3)

	t = [cost.sps cost.sd cost.pcg cost.os8 cost.pre];
%	disp(t(1:20,:))
end
