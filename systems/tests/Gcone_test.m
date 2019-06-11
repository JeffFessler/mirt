% Gcone_test.m
% test the Gcone object
% Copyright 2008-1-1, Jeff Fessler, University of Michigan

printm 'todo: check user_source_zs vs hct2 with file_source_z'

pn = jf_protected_names;

f.sf1 = {'sf1', ...
	'sf1,p:ts,b:st', 'sf1,p:st,b:st', ...
	'sf1,p:ts,b:ts', 'sf1,p:st,b:ts'};
f.sf2 = {'sf2', ...
	'sf2,p:ts,b:st', 'sf2,p:st,b:st', ...
	'sf2,p:ts,b:ts', 'sf2,p:st,b:ts'};
systypes = { ... % lots of variations on system models!
	'rf1', ... % RR model
	f.sf1{:}, f.sf2{:}, ... % SF variations
	'sf1', 'sf2', 'sf3', 'sf4', 'sf5' ...
};
%systypes = {'rf1', 'sf1', 'dd2'}
if 0 && exist('dd_ge1_mex') == 3 % UM only
	systypes{end+1} = 'dd1';
	systypes{end+1} = 'dd2'; % todo
end
% systypes = {'nn1', 'pd1'};
%systypes = {'dd2'}; % todo! fails adjoint test!?
%systypes = {'sf1', 'sf2'};
%systypes = f.sf2;
nsys = numel(systypes);

%f.class = 'Fatrix';
f.class = 'fatrix2';

% todo: test both ways for fatrix2, not for Fatrix
%is_zxy = false;
is_zxy = true;
if is_zxy
	to_zxy = @(x) permute(x, [3 1 2]);
else
	to_zxy = @(x) x;
end

% small systems for basic tests
if 1 && (0 || ~isvar('A1')), printm 'setup small'
	f.down = 16;

	dfs_list = [inf inf 0]; % parallel flat arc
	dsd_list = [inf 949.075 949.075]; % parallel flat arc

	if 0 % arc only
		dfs_list = [0];
		dsd_list = [949.075];
	end

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

			if streq(systype, 'dd1') || streq(systype, 'dd2')
				f.dy = 1; % DD requires square pixels
				if isinf(cgs.dsd) % no parallel for DD
					continue
				end
			else
				f.dy = 0.7; % stress test SF with non-square
			%	f.dy = 'dx';
			end

			igs = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
				'dx', 1, 'dy', f.dy, 'dz', 0.5, ...
				'offset_x', 2.9, 'offset_y', 3.5, ...
				'offset_z', -3.4, ...
				'mask', 'all-but-edge', ... % trick: for hct2
				'down', f.down);

			if im, cgs.plot(igs); end

			args = {cgs, igs, 'type', systype, 'class', f.class};
			if streq(systype, 'sf4') || streq(systype, 'sf5')
				args = {args{:}, 'mexarg', {single(0)}};
			end

			args = {args{:}, 'zxy', is_zxy};

			A1 = Gcone(args{:}, 'nthread', 1);
			Ac = Gcone(args{:});

			if 1 % adjoint
				tester_tomo2(A1, to_zxy(igs.mask), 'G2', Ac)
				test_adjoint(A1, 'big', 1, 'tol', 5e-5)
				test_adjoint(Ac, 'big', 1, 'tol', 5e-5)
			end

			if streq(systype, 'dd1') || streq(systype, 'dd2')
				continue % don't bother hct2 test for DD
			end

			% hereafter is hct2 test
			Ah = Gcone(args{:}, 'use_hct2', 1);

			if pn.has_hct2
				if streq(systype, 'nn1') || streq(systype, 'pd1')
					thresh = 3e-2; % big due to rounding in nn1
				else
					thresh = 4e-6;
				end

				if 1 % mex vs hct2
					xs = single(igs.mask);
					xs(round(end/3), round(2*end/3), round(end/4)) = 10;
					xs = to_zxy(xs);
					t1 = A1 * xs;
					t2 = Ah * xs;
					equivs(t1, t2, 'thresh', thresh)

					b1 = A1' * t1;
					b2 = Ah' * t1;
%					im plc 1 3, im(1, b1), im(2, b2), im(3, b2-b1)
					equivs(b1, b2, 'thresh', thresh)
				end
				if 0 % block
					B1 = Gblock(A1, 2);
					Bh = Gblock(Ah, 2);
					t1 = B1{2} * xs;
					t2 = Bh{2} * xs;
					equivs(t1, t2, 'thresh', thresh)
				end

				tester_tomo2(Ah, to_zxy(igs.mask), ...
					'equiv_thresh', thresh) % paces
				test_adjoint(Ah, 'big', 1, 'tol', 5e-5)
			end
		end % systype
	end % dfs
end % small


if 0, printm 'rf1' % test rf1 vs dd2 small
	A0 = Gcone(cgs, igs, 'type', 'rf1');
	Ad = Gcone(cgs, igs, 'type', 'dd2');
%	x0 = single(igs.mask);
	x0 = igs.circ;
	im(x0)
	cpu etic
	yd = Ad * x0;
	cpu etoc dd
	if 1
		cpu etic
		y0 = A0 * x0;
		cpu etoc rf1
	end
	im plc 2 2
	im(yd), im(y0), im(y0-yd)
return
end


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
%	clf, im(x0), return
%prompt
end


% big systems for accuracy tests
if 0 || ~isvar('Ab'), printm 'setup big'
	cgb = ct_geom('ge1', 'nt', 320, ...
		'source_z0', -40, 'pitch', 0.5, ... % test helix
	...%	'dfs', inf, ... % flat detector
	...%	'dsd', inf, 'dfs', inf, ... % parallel beam
		'down', f.down);

	clear Ab Ah
	for ii=1:nsys
		systype = systypes{ii};
		if streq(systype, 'dd', 2) && isinf(cgb.dsd)
			Ab{ii} = [];
			continue
		end
		Ab{ii} = Gcone(cgb, igb, 'type', systype, ...
				'zxy', is_zxy, 'class', f.class);
	end
end


if ~isvar('ya'), printm 'analytical projections'
	ya = ellipsoid_proj(cgb, ell, 'oversample', 2);
%	im clf, im(ya)
%prompt
end


if 0 || ~isvar('yb'), printm 'discrete projections'
	nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
	for ii=1:nsys
		if ~isempty(Ab{ii})
			cpu etic
			yb{ii} = Ab{ii} * to_zxy(x0);
			f.time(ii) = cpu('etoc');
			printm('%s: time %g nrmse %g %%', ...
				systypes{ii}, f.time(ii), nrmse(ya, yb{ii}))
		else
			f.time(ii) = 0;
			yb{ii} = cgb.zeros;
		end
	end
end


if 0 % toggle analytical vs big discrete
	ii = 1;
	im_toggle(ya(:,end/2,:), yb{ii}(:,end/2,:), ...
		yb{ii}(:,end/2,:) - ya(:,end/2,:))
	im(ya(:,end/2,:) - yb{ii}(:,end/2,:))
return
end


if 0, printm 'look at error in worst views'
	ilist = 1:nsys;
	ilist = [1 2 7];
	im('plc', 2, numel(ilist))
	for jj=1:numel(ilist)
		ii = ilist(jj);
		err = yb{ii} - ya;
		tmp = reshape(err, [], cgb.na);
		tmp = sum(tmp.^2); % error in each view
		ia = imax(tmp); % worst view
		im(jj, err(:,:,ia)), cbar h
		titlef('%s ia=%d', systypes{ii}, ia)
		im('subplot', jj+numel(ilist))
		plot(tmp), axis tight
	end
return
end


if 1, printm 'projection profiles'
	it = cgb.nt;
	it = round(cgb.nt/2); % middle
	it = it + [-2 0 2];
	ia = imin(abs(cgb.ad - 45)); % closest to 45
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
		ir_legend({systypes{:}, 'true'})
		axisy(0, 1.2 * max(col(pro(ya))))
		grid
	end
%prompt
return
end


if 0 % dd1 vs dd2 - they match well
	i_dd1 = strmatch('dd1', systypes);
	i_dd2 = strmatch('dd2', systypes);
	im clf, im(yb{i_dd1} - yb{i_dd2}), cbar
	equivs(yb{i_dd1}, yb{i_dd2}, 'thresh', 2e-5)
return
end


if 1, printm 'show projections and differences'
	im('plc', 2, 1+nsys)
	im(1, x0)
	ia = round([1 cgb.na/4 cgb.na/2 cgb.na]);
	im(nsys+2, ya(:,:,ia)), axis normal

	for ii=1:nsys
		tmp = yb{ii};
		im(ii+1, tmp(:,:,ia)), axis normal
		xlabel(systypes{ii})
		im(ii+2+nsys, tmp(:,:,ia)-ya(:,:,ia)), axis normal
	end
prompt
end


if 1, printm 'show back-projections and pairwise differences'
	im('plc', nsys, nsys)
	iz = round([1 igb.nz/4 igb.nz/2 igb.nz]);
	iz = 1:2:igb.nz;
	clear yt
	for ii=1:nsys
%		tmp = cgb.ones;
		tmp = cgb.zeros; tmp(:,:,20) = 1;
		yt{ii} = Ab{ii}' * tmp;
	end
	for ii=1:nsys
		for jj=1:nsys
			if ii == jj
				im('subplot', (ii-1)*nsys + jj)
				im(yt{ii}(:,:,iz))
				xlabel(systypes{ii})
			elseif ii > jj
				im('subplot', (ii-1)*nsys + jj)
				im(yt{ii}(:,:,iz) - yt{jj}(:,:,iz))
				xlabelf([systypes{ii} ' - ' systypes{jj}])
			end
		end
	end
end
