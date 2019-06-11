% Reg1_test.m
% test of Reg1 regularization object

if 1 % test single-pixel case
	R0 = Reg1(true(1));
	x = 7;
	jf_equal(0, R0.penal(R0, x))
	jf_equal(0, R0.cgrad(R0, x))
	jf_equal(0, R0.denom(R0, x))
	jf_equal(full(R0.C), nan(0,1))
end


if 1 % test tiny image
	R0 = Reg1(true(3,2), 'offsets', '2d:hv');
	x = reshape(1:3*2, 3, 2);
	R0.penal(R0, x);
	R0.cgrad(R0, x);
	R0.denom(R0, x);
	full(R0.C);

	tmp = {true(3,2), 'offsets', '2d:hv'};
	if 0
		rng(7)
		tmp = {tmp{:}, 'user_wt', ones(3,2,2)};
	end
	R0t = Reg1(tmp{:}, 'type_penal', 'mat', 'type_diff', 'ind');
	R0x = Reg1(tmp{:}, 'type_penal', 'mex');
	jf_equal( full(R0t.C), full(R0x.C) )
	jf_equal( full(R0t.C1), full(R0x.C1) )
%	jf_equal( R0t.cgrad(R0t, x), R0x.cgrad(R0x, x) ) % not ?
%	full(R0t.C1)
%	full(R0x.C1)
end


if 1 || ~isvar('Rt'), printm 'Rt'
	if 1 % 3d
		ig = image_geom('nx', 512, 'ny', 480, 'nz', 2^5, 'fov', 500, ...
			'down', 2^2);
	else % 2d
		ig = image_geom('nx', 512, 'ny', 448, 'fov', 500, 'down', 2^3);
	end
	tmp = ig.circ(ig.fov/2*1.1) > 0;
	tmp([1 end],:,:) = 0; tmp(:, [1 end],:) = 0; % zero border for 3d
	ig.mask = tmp; clear tmp
	kappa = 2 + 1 * cos(2*pi*ig.x/ig.fov) * sin(2*pi*ig.y/ig.fov)';
%	kappa = ig.mask;
	f.offsets = ''; f.n_offset = 4;
	if ig.is3
		kappa = repmat(ig.mask_or, [1 1 ig.nz]) .* ig.mask;
		f.offsets = '3d:26';
		f.n_offset = 13;
	else
		kappa = single(kappa) .* ig.mask;
	end

	if 0
		profile on
		cpu etic
		pn = jf_protected_names;
		for qq=1:29*2
			ii = 1:prod(ig.dim);
%			[jx jy jz] = ind2sub(ig.dim, ii);
			jj = pn.ind2sub(ig.dim, ii);
		end
		cpu etoc
		profile report
	return
	end

	f.l2b = 3;
%	f.pot = 'hyper3'; f.delta = 1; f.pot_arg = {f.pot, f.delta};
%	f.pot = 'quad'; f.delta = inf; f.pot_arg = {f.pot};
%	f.pot = 'qgg2'; f.delta = 2; f.pot_arg = {f.pot, f.delta, 1.2};
	f.pot = 'gf1'; f.delta = 2; ...
		tmp = potential_fun('gf1-fit', nan, {'qgg2', [1 10], 1.2});
		f.pot_arg = {f.pot, f.delta, tmp};

%	x = ig.unitv(ig.nx, ig.ny) + ig.unitv(ig.nx/4, 5+0*ig.ny/2);
	x = single(ig.mask + ig.circ(ig.fov/4));
	rng(0)
	x = single(rand(ig.dim)) .* ig.mask;
	if ig.is3
		x = x + ig.unitv ...
			+ ig.unitv(3*ig.nx/8, ig.ny/2, 1) ...
			+ ig.unitv(ig.nx/2, 5*ig.ny/8, ig.nz);
	end
	xm = x(ig.mask);

	for order = 1:2
		pr order
		f.arg1 = {kappa, 'offsets', f.offsets, 'beta', 2^f.l2b, ...
			'edge_type', 'tight', 'order', order};
		if 0 && order == 1 % test user_wt
			rng(7)
%			tmp = ones([ig.dim f.n_offset]);
%			tmp(1:end/4,:,:) = 4;
			tmp = rand([ig.dim f.n_offset]);
			f.arg1 = {f.arg1{:}, 'user_wt', tmp};
		end
		f.arg = {f.arg1{:}, 'pot_arg', f.pot_arg};
		Rt = Reg1(f.arg{:}, 'type_penal', 'mat');
		Rx = Reg1(f.arg{:}, 'type_penal', 'mex'); % default control
%		R1 = Reg1(f.arg{:}, 'type_penal', 'mex', 'control', 1);

		if 0 % test timing of loop control 1 and 2
			% in 3d for order=2 with no user_wt,
			% control=2 is about 5% faster
			g1 = ig.embed(R1.cgrad(R1, xm)); % warm up
			g2 = ig.embed(Rx.cgrad(Rx, xm)); % warm up

			cpu etic
			g1 = ig.embed(R1.cgrad(R1, xm));
			cpu etoc 'Rx cgrad control=1 time'

			cpu etic
			g2 = ig.embed(Rx.cgrad(Rx, xm));
			cpu etoc 'Rx cgrad control=2 time'
			equivs(g1, g2)
		continue
		end

		if has_mex_jf
			switch f.pot
			case {'qgg2', 'gf1'}
				tmp = {f.delta, f.pot_arg{3}};
			otherwise
				tmp = f.delta;
			end
			Ro = Robject(f.arg1{:}, 'type_denom', 'matlab', ...
				'potential', f.pot, 'delta', tmp, ...
				'distance_power', Rt.distance_power);
		end

		if 1 && ig.is3 && has_mex_jf
			zxy = @(x) permute(x, [3 1 2]);
			xyz = @(x) permute(x, [2 3 1]);
 			% convert offsets to zxy!
			f.offsets_zxy = reg_offset_xyz_to_zxy(f.offsets, size(kappa));
			Rz = Reg1(zxy(kappa), f.arg{2:end}, ...
				'user_wt', [], ...
				'offsets', f.offsets_zxy, ...
				'offsets_is_zxy', true, ...
				'type_penal', 'zxy');
		end

%		Rm = Rmem(kappa, 'pot_arg', f.pot_arg, 'beta', 2^f.l2b, ...
%			'edge_type', 'tight', ...
%			'distance_power', Rg.distance_power, ...
%			'offsets', f.offsets);

		if 1 % check C1 *
			jf_equal(Rt.C1 * x, Rx.C1 * x)
			jf_equal(Rt.C1 * xm, Rx.C1 * xm)
			jf_equal(ig.shape(Rt.C1 * xm), Rt.C1 * x)
			if has_mex_jf
				jf_equal(Ro.C1 * x, Rx.C1 * x)
			end
		end

		if 0 % check C1 (only for tiny cases)
			c1 = Rt.C1; c1 = c1(:,:);
			c2 = Rx.C1; c2 = c1(:,:);
			jf_equal(c1, c2)
			if has_mex_jf
				co = Ro.C1; co = co(:,:);
				jf_equal(c1, co)
			end
%			clear c1 c2
		end

		if 1 % check dercurv
			[t.dt t.ct] = feval(Rt.dercurv, Rt, Rt.C1 * x);
			[t.dx t.cx] = feval(Rx.dercurv, Rx, Rx.C1 * x);
%			[t.d1 t.c1] = feval(R1.dercurv, R1, R1.C1 * x);
			jf_equal(t.dt, t.dx)
			jf_equal(t.ct, t.cx)
%			jf_equal(t.dt, t.d1)
%			jf_equal(t.ct, t.c1)

			[t.dx t.cx] = feval(Rx.dercurv, Rx, Rx.C1 * xm);
			jf_equal(ig.shape(t.dx), t.dt)
			jf_equal(ig.shape(t.cx), t.ct)

			if has_mex_jf
				[t.do t.co] = feval(Ro.dercurv, Ro, Ro.C1 * xm);
				t.do = ig.shape(t.do);
				t.co = ig.shape(t.co);
				equivs(t.dt, t.do)
				% trick: curv matches *except* at outer edge(s)
				equivs(t.ct(2:end-1,2:end-1,:), t.co(2:end-1,2:end-1,:))
			end
		end

		if 1 && streq(f.pot, 'quad') % check C *
%			tmp = Rt.C * x;
			jf_equal(ig.shape(Rt.C * xm), Rt.C * x)
		end

		if 1 % check penalty value
			tmp = Rt.penal(Rt, x);
			jf_equal(tmp, Rt.penal(Rt, xm))
			equivs(tmp, Rx.penal(Rx, x))
%			equivs(tmp, R1.penal(R1, x))
			jf_equal(Rx.penal(Rx, x), Rx.penal(Rx, xm))
			if has_mex_jf
%				jf_equal(Rt.penal(Rt, xm), Ro.penal(Ro, xm))
				equivs(tmp, Ro.penal(Ro, xm), 'thresh', 18e-6)
			end
		end

		if 1 % check cgrad
			g1 = ig.embed(Rt.cgrad(Rt, xm));
			w1 = ig.embed(Rt.denom(Rt, xm));
			w1s = ig.embed(Rt.denom_sqs1(Rt, xm));
			if 1
				cpu etic
				g2 = Rt.cgrad(Rt, x);
				cpu etoc 'Rt cgrad time'
				jf_equal(g1, g2)
			end

			if has_mex_jf
				go = ig.embed(Ro.cgrad(Ro, xm));
				equivs(g1, go)

				g3 = ig.embed(Rx.cgrad(Rx, xm));
				equivs(g1, g3)

				if 1 % check mex version w/o mask
					cpu etic
					g4 = Rx.cgrad(Rx, x);
					cpu etoc 'Rx cgrad time'
					jf_equal(g3, g4)
				end

				if ig.is3 && order == 1
					tmp = zxy(x);
					% check zxy version w/o mask
					cpu etic
					[g5 w5] = feval(Rz.cgrad_denom, Rz, tmp);
					cpu etoc 'Rz cgrad/denom time'
					g5 = xyz(g5);
					w5 = xyz(w5);
					if isempty(Rx.user_wt) % trick
						equivs(g1, g5)
						equivs(w1, w5)
					end
%					im(w1-w5), return

					% zxy version with mask
					mask_zxy = zxy(ig.mask);
					[g5m w5m] = feval(Rz.cgrad_denom, Rz, tmp(mask_zxy));
					g5m = xyz(embed(g5m, Rz.mask));
					w5m = xyz(embed(w5m, Rz.mask));
					jf_equal(g5, g5m)
					jf_equal(w5, w5m)
				end
			end

			if 0 % check empirical gradient (slow)
				ge = ig.embed(Rx.egrad(Rx, xm, 0.001));
				max_percent_diff(ge, g1)
				im plc 2 3
				im(1, x, 'x'), cbar
				im(2, g1, 'g mat'), cbar
				im(3, g3, 'g mex'), cbar
				im(4, ge, 'g emp'), cbar
				im(5, g1-ge, 'mat-emp'), cbar
				im(6, g3-ge, 'mex-emp'), cbar
			prompt
			end

%			im(g3(:,:,[1 end])-g1(:,:,[1 end]))
%			im clf, im([g1; g3; g1-g3]), cbar

%			clear g1 g2 g3 g4
		end

		if 1 % check denom
			dentf = Rt.denom_sqs1(Rt, x);
			dentm = Rt.denom_sqs1(Rt, xm);
			dentm = ig.embed(dentm);
			equivs(dentf, dentm)

			if has_mex_jf
				denxf = Rx.denom_sqs1(Rx, x);
				denxm = Rx.denom_sqs1(Rx, xm);
				denxm = ig.embed(denxm);
				equivs(dentf, denxf)
				equivs(denxm, denxf)
%				im clf, im([denxf; denxm; denxf-denxm]), cbar

%				denx1 = R1.denom_sqs1(R1, x);
%				equivs(dentf, denx1)

				dentf = Rt.denom(Rt, x);
				denxf = Rx.denom(Rx, x);
%				denx1 = R1.denom(R1, x);
				equivs(dentf, denxf)
%				equivs(dentf, denx1)

				deno = ig.embed(Ro.denom(Ro, xm));
				% trick: denom matches except along outer edges
				ix = (1+order):(ig.nx-order);
				iy = (1+order):(ig.ny-order);
				equivs(dentf(ix,iy), deno(ix,iy))

				if order == 1 && ig.is3 && isempty(Rx.user_wt)
					equivs(dentf, w5) % check zxy
				end
			end
		end

		if 0 % check diag ?
			%d = Rm.diag(R);
		end

		% test threads
	end
prompt
end


return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robject designed to match Rmem'
if ~isvar('Ro'), printm 'Ro'
	tmp = repmat(kappa.^2, [1 1 1 length(Rm.offsets)]);
	Ro = Robject(ig.mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
		'offsets', Rm.offsets, ...
		'potential', f.pot, 'delta', f.delta, ...
		'type_denom', 'matlab', ...
		'distance_power', 0, 'user_wt', tmp);
	clear tmp
end

cpu etic
g1 = Rm.cgrad(Rm, x);
cpu etoc 'Rm cgrad'
cpu etic
g2 = Ro.cgrad(Ro, x(ig.mask));
g2 = embed(g2, ig.mask);
cpu etoc 'Ro cgrad'
max_percent_diff(g1, g2)
if ~isequal(g1, g2), error 'cgrad mismatch', end, clear g2
