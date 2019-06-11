% Gcone_twv_test.m
% test the Gcone object for "thread within view" version by Brian Gonzalez
% Copyright 2015-07-20, Jeff Fessler, University of Michigan

pn = jf_protected_names;

systypes = 'sf2,v:all';
systypes = {'sf2', systypes}
nsys = numel(systypes);

f.is_zxy = true;
if f.is_zxy
	to_zxy = @(x) permute(x, [3 1 2]); to_xyz = @(z) permute(z, [2 3 1]);
else
	to_zxy = @(x) x;
end

% small systems for basic tests
if 1 && (0 || ~isvar('Ac')), printm 'setup small'
	f.down = 16;

	dfs_list = [0 inf inf]; % arc flat parallel
	dsd_list = [949.075 949.075 inf]; % arc flat parallel

	f.dy = 0.7; % stress test SF with non-square
	igs = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dy', f.dy, 'dz', 0.5, ...
		'offset_x', 2.9, 'offset_y', 3.5, ...
		'offset_z', -3.4, ...
		'mask', 'all-but-edge', ... % trick: for hct2
		'down', f.down);

	for kk = 1:numel(dfs_list)

		cgs = ct_geom('ge2', 'nt', 320, ...
			'source_z0', -20, 'pitch', 0.5, ... % test helix
			'dfs', dfs_list(kk), ... % arc or flat
			'dsd', dsd_list(kk), ... % fan or parallel beam
			'down', f.down);

		for ii=1:nsys
			systype = systypes{ii};
			printm('testing type %s dfs=%g dsd=%g', ...
				systype, cgs.dfs, cgs.dsd)

			args = {cgs, igs, 'type', systype};
			args = {args{:}, 'zxy', f.is_zxy};

			Ac = Gcone(args{:});

		%	im(Ac * to_zxy(single(igs.mask)))
		%	im(to_xyz(Ac' * cgs.ones)), prompt

			tester_tomo2(Ac, to_zxy(igs.mask)) % paces
			test_adjoint(Ac, 'big', 1, 'tol', 5e-5)
		end % systype
	end % dfs
end % small


if ~isvar('x0'), printm 'x0 big'
	f.down = 4;
	igb = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
		'dy', 0.7, ... % stress test
		'offset_x', 12.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'down', f.down);
	ell = [3*igb.dx 5*igb.dx -2*igb.dz ...
		igb.dx*igb.nx/3 igb.dy*igb.ny/4 igb.zfov/4 ...
		0 0 10];
	x0 = ellipsoid_im(igb, ell, 'oversample', 2);
end


% big systems for accuracy tests
if 0 || ~isvar('Ab'), printm 'setup big'
	cgb = ct_geom('ge1', 'nt', 320, ...
		'source_z0', -40, 'pitch', 0.5, ... % test helix
		'down', f.down);
%		'dfs', inf, ... % flat detector
%		'dsd', inf, 'dfs', inf, ... % parallel beam

	for ii=1:nsys
		systype = systypes{ii};
		Ab{ii} = Gcone(cgb, igb, 'type', systype, 'zxy', f.is_zxy);
	end
end


if ~isvar('ya'), printm 'analytical projections'
	ya = ellipsoid_proj(cgb, ell, 'oversample', 2);
%	im clf, im(ya)
%prompt
end


if 0 || ~isvar('yb'), printm 'discrete projections'
	nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
	x00 = to_zxy(x0);
	for ii=1:nsys
		cpu etic
		yb{ii} = Ab{ii} * x00;
		f.time(ii) = cpu('etoc');
		printm('%s: time %g nrmse %g %%', ...
			systypes{ii}, f.time(ii), nrmse(ya, yb{ii}))
	end
end


if 0
	im_toggle(ya(:,end/2,:), yb{1}(:,end/2,:), ...
		yb{1}(:,end/2,:) - ya(:,end/2,:))
	im(ya(:,end/2,:) - yb{1}(:,end/2,:))
	im(yb{1}(:,end/2,:) - yb{2}(:,end/2,:))
end

if 1, printm 'look at error in worst views'
	im('pl', 2, nsys+1)
	for ii=1:numel(systypes)
		err = yb{ii} - ya;
		tmp = reshape(err, [], cgb.na);
		tmp = sqrt(mean(tmp.^2)); % error in each view
		ia = imax(tmp); % worst view
		im(ii, err(:,:,ia)), cbar h
		titlef('%s ia=%d', systypes{ii}, ia)
		im('subplot', ii+1+nsys)
		if im
			plot(1:cgb.na, tmp), axis tight
			titlef('%s RMS', systypes{ii})
			xlabelf('view')
		end
	end
	if 1
		err = yb{2} - yb{1};
		tmp = reshape(err, [], cgb.na);
		tmp = sqrt(mean(tmp.^2)); % error in each view
		ia = imax(tmp); % worst view
		im(nsys+1, err(:,:,ia)), cbar h
		titlef('%s vs %s', systypes{1}, systypes{2})
	end
return
end

if 1, printm 'projection profiles'
	it = cgb.nt;
	it = round(cgb.nt/2); % middle
	it = it + [-2 0 2];
	ia = imin(abs(cgb.ad - 45));
%	ia = ceil(ia/2);
	pro = @(y) col(y(:,it,ia));
	arg = [];
	for ii=1:numel(systypes)
		arg = [arg pro(yb{ii})];
	end
	if im
		clf, plot([arg pro(ya)])
		text(10, 200, sprintf('ia=%d', ia))
		text(10, 400, sprintf('ang=%g', cgb.ad(ia)))
		legend(systypes{:}, 'true')
		axisy(0, 1.2 * max(ya(:)))
		grid
	end
return
end
