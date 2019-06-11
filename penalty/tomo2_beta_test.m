% tomo2_beta_test
% Examine how the beta that is required to ensure a certain resolution
% changes as a function of tomographic sampling parameters.
% Summary: need to multiply beta by na/dr to keep consistent resolution.
% For somewhat consistent resolution in *pixels* we would also scale beta
% by dx^3, but we want consistent resolution in *mm* and no scaling w.r.t.
% dx is need for that, due to rho^3 term in local Fourier analysis.
% Copyright 2005-8-26, Jeff Fessler, The University of Michigan

redo = 0;
if redo || ~isvar('G2'), disp 'setup G'
	down = 4;
	ig1 = image_geom('nx', 512, 'ny', 272, 'fov', 256, 'down', down);
	sg1 = sino_geom('par', 'nb', 888, 'na', 984, 'dr', 541/949, ...
		'down', down);
	ig2 = ig1;
	sg2 = sg1;

	% vary the next three to investigate their effect
	sg2.na = sg2.na / 1;
	sg2.dr = sg2.dr / 1;
	dscale = 2;
	ig2.dx = ig2.dx / dscale;
	ig2.dy = ig2.dy / dscale; ig2.fov = ig2.fov / dscale;

	f.mask = [test_dir 'mask.fld'];
	fld_write(f.mask, ig1.mask, 'type', 'byte');

	if 0 % linear interpolation (cf system 9)
		f.system = 9;
		table = {'linear', 'chat', 0, 'Ntab', 1000};

	elseif 1 % square-pixel / strip-integral
		f.system = 2;

		if 0
			f.strip_width = 0;
		else
			f.strip_width = sg1.dr;
		end
		f.table = {'square/strip', 'chat', 0, 'Ltab', 1000, ...
			'strip_width', f.strip_width};
	end

	if 1
		G1 = Gtomo2_table(sg1, ig1, f.table, 'nthread', 2);
		G2 = Gtomo2_table(sg2, ig2, f.table, 'nthread', 2);

		if 0, printm 'table vs strip'
			G2 = Gtomo2_strip(sg1, ig1, 'strip_width', f.strip_width);
			im pl 1 3
			y1 = G1 * ig1.mask;
			y2 = G2 * ig1.mask;
			y1 = G1' * sg1.unitv(60, 30);
			y2 = G2' * sg1.unitv(60, 30);
			im(1, y1), cbar
			im(2, y2), cbar
			im(3, y2-y1), cbar
		return
		end
	else
		f.arg1 = aspire_pair(sg1, ig1, 'strip_width', f.strip_width, ...
			'support', ['file ' f.mask], 'system', f.system);
		f.arg2 = aspire_pair(sg2, ig2, 'strip_width', f.strip_width, ...
			'support', ['file ' f.mask], 'system', f.system);
		G1 = Gtomo2_dscmex(f.arg1);
		G2 = Gtomo2_dscmex(f.arg2);

%	else
%		G1 = Gtomo2_dd(mask, f.arg{:});
	end
prompt
end

if 0 % check scales of A'A
	f.offset = round([ig1.nx ig1.ny] .* [.07 .1]);
	ej = ig1.unitv(ig1.nx/2+1+f.offset(1),ig1.ny/2+1+f.offset(2));
	aa1 = G1' * (G1 * ej(ig1.mask));
	aa2 = G2' * (G2 * ej(ig1.mask));
	aa1 = aa1 / sg1.na / ig1.dx^3 * sg1.dr;
	aa2 = aa2 / sg2.na / ig2.dx^3 * sg2.dr;
	sum([aa1(:) aa2(:)])
	t = [aa1 aa2];
	im clf, im(t), cbar
return
end

if redo || ~isvar('mtf1')
	beta0 = 2^3;
	R0 = Robject(ig1.mask, 'beta', beta0);
%	R0 = Reg1(ig1.mask, 'beta', beta0);
	scale1 = sg1.na * ig1.dx^0 / sg1.dr; % this is the key scaling
	scale2 = sg2.na * ig2.dx^0 / sg2.dr;
	o1 = 0*round([6/ig1.dx 8/ig1.dy]);
	o2 = 0*round([6/ig2.dx 8/ig2.dy]);
	[psf1 var1 fw1 mtf1] = qpwls_psf(G1, R0, scale1, ...
		ig1.mask, 1., 'offset', o1, 'dx', ig1.dx);
	[psf2 var2 fw2 mtf2] = qpwls_psf(G2, R0, scale2, ...
		ig2.mask, 1., 'offset', o2, 'dx', ig2.dx);
	printf('var = %g %g', var1, var2)
	if im
		im clf
		im(121, ig1.x, ig1.y, psf1), cbar
		axis equal
		axis([6+[-1 1]*9 8+[-1 1]*9])
		im(122, ig2.x, ig2.y, psf2), cbar
%		axis([-1 1 -1 1]*10)
		axis equal
		axis([6+[-1 1]*9 8+[-1 1]*9])
	%	im(212, psf1-psf2), cbar
	end
end

if 1 && im % plot profiles
	clf, subplot(211)
	ix = ig1.dx * (o1(1) + 0.5);
	iy = o1(2) + ig1.ny/2+1;
	plot((ig1.x - ix)/4, psf1(:,iy) / max(psf1(:)), '-oy', 0, 1, '.-c')
	axis([[-12 12]/4 -0.05 1.05])
	legend('\Delta_x = 1', '\Delta_x = 1/2')
	title('Profiles through PSFs for same \alpha')
	xlabel 'spatial position [mm]'

	ix = ig2.dx * (o2(1) + 0.5);
	iy = o2(2) + ig2.ny/2+1;
	hold on
	plot((ig2.x - ix)/4, psf2(:,iy) / max(psf2(:)), '.-c')
	hold off
%	ir_savefig c psf_profile
end
