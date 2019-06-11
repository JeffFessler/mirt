% tml_os_mltr_test.m
% test ML-TR (OS) with different preconditioners (denominators)
% Copyright 2010-07-11, Jeff Fessler, University of Michigan

% generate data
if ~isvar('xfbp'), printm 'setup tml_os_mltr_test'
	ig = image_geom('nx', 120, 'fov', 25 * 0.2 / 1000); % trick: HU!
	ig.mask = ig.circ > 0;
	[xtrue ell] = ellipse_im(ig, [], 'oversample', 2, 'hu_scale', 1000);
	im plc 3 2
	f.clim = 1000 + [-1 1] * 50;
	im(1, xtrue, f.clim), cbar

	if 1 % parallel-beam
		sg = sino_geom('par', 'nb', 120, 'na', 1*120, ...
			'dr', ig.dx, 'strip_width', 'd');
		A = Gtomo2_wtmex(sg, ig, 'nthread', jf('ncore'), 'chat', 0);
	else
		sg = sino_geom('moj', 'nb', 240, 'na', 120, ...
			'dx', ig.dx, 'strip_width', 'd');
		A = Gtomo2_table(sg, ig, {'mojette,square/strip'}, ...
			'nthread', jf('ncore'))
	end

%	proj = A * xtrue;
	proj = ellipse_sino(sg, ell, 'oversample', 3);
%	minmax(proj)

	bi = 10 * sg.ones;
	ytrue = bi .* exp(-proj);
	ri = sg.zeros;
	yi = ytrue; % noiseless for now

	im(3, ytrue), cbar
	im(4, yi), cbar

%	xfbp = (1/f.scale) * tr_fbp(sg, ig, yi, bi, ri);
	xfbp = tr_fbp(sg, ig, yi, bi, ri);
	im(2, xfbp, f.clim), cbar
prompt
end

if ~isvar('Ab'), printm 'Ab'
	f.nblock = 40;
	Ab = Gblock(A, f.nblock);
end

% run ML-TR
if ~isvar('xmltr'), printm 'matlab T-ML-OS-SPS'
	f.niter = 11;
	f.step = 1;
%	xinit = ig.ones * mean(xfbp(ig.mask));
	xinit = ig.ones * 1000; % water

	% try each precon type
	plist = {'classic', 'fast'}; % 'nuyts'

	for ii = 1:length(plist)
		tmp = tml_os_mltr(xinit(ig.mask), Ab, yi, bi, ri, ...
			'precon', plist{ii}, ...
			'step', f.step, ...
			'niter', f.niter, 'isave', 'all');
		xmltr{ii} = ig.embed(tmp);
		im clf, im(xmltr{ii}, f.clim), cbar
		titlef('T-ML-OS-TR iterations: %s', plist{ii})
	prompt
	end

	max_percent_diff(xmltr{1}, xmltr{2})
end

if 1
	im plc 3 2
	im(1, xtrue, 'True', f.clim), cbar 
	im(2, xfbp, 'FBP', f.clim), cbar 

	for ii = 1:length(plist)
		im(2+ii, xmltr{ii}(:,:,end), f.clim), cbar
		titlef('T-ML-OS-TR iterations: %s', plist{ii})
		im(4+ii, xmltr{ii}, f.clim), cbar
%		im(4+ii, xmltr{ii}), cbar
	end
end
