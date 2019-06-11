% feldkamp_example.m
% Example of how to use feldkamp.m for cone-beam CT reconstruction
% Copyright 2004-8-28, Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, University of Michigan

if ~isvar('cg'), printm 'cg: cone-beam CT geometry'
	% see book chapter section (ask me) on cone-beam CT recon for notation
	if ~isvar('down'), down = 4; end % down sample a lot to save time
	% default is arc detector; but allow flat for cone_beam_ct_example.m
	if ~isvar('dfs'), dfs = 0; end
	cg = ct_geom('fan', 'ns', 256, 'nt', 240, 'na', 288, ...
		'ds', 1024/256, 'dt', 1024/256, ...
		'down', down, ...
		'offset_s', 0.25, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', 949, 'dod', 408, 'dfs', dfs);
	printm('fov rmax=%g', cg.rmax)
	clear dfs
end


if ~isvar('ig'), printm 'ig: image geometry'
	ig = image_geom('nx', 256, 'ny', 240, 'nz', 200, 'fov', 500, ...
		'down', down);
%	       'dz', -6); % negative dz to match aspire
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2
end


if ~isvar('ell'), printm 'ell: ellipsoid object'
	cyl = [20 10 10  150 150 inf  0 0.01]; % 30cm diam cylinder
	cyl(end) = 0; % todo: kludge until cylinder_proj is fixed
	ell = [ ...
		[20 10 10	150 150 980	0 0 0.01]; % 30cm diam "cylinder
		[80 10 10	50 50 30	0 0 0.01]; % bone-like inserts
		[-10 -40 75	40 40 40	0 0 0.01];
		[-10 80 -20	30 30 30	0 0 0.01];
	];
end


if ~isvar('xtrue'), printm 'xtrue: true image volume'
	f.over = 2;
	tmp = ellipse_im(ig, cyl(:,[1 2 4 5 7 8]), 'oversample', f.over);
	tmp = repmat(tmp, [1 1 ig.nz]); % emulate a "cylinder_im"
	xtrue = tmp + ellipsoid_im(ig, ell, 'oversample', f.over);

	im plc 2 3
	im(1, ig.x, ig.y, xtrue), cbar
	titlef('x true, z=%g to %g', ig.z(1), ig.z(end))
	clear tmp
end


if ~isvar('proj'), printm 'proj: analytical ellipsoid projection views'
	tmp = cylinder_proj(cg, cyl, 'oversample', f.over);
	proj = tmp + ellipsoid_proj(cg, ell, 'oversample', f.over);
	im(4, proj, 'true projections'), cbar
	clear tmp
prompt
end


% noisy data and estimated line integrals
if ~isvar('li_hat'), printm 'li_hat: projection view log data'
	% noisy data, if blank scan value has been specified.
	if isvar('bi') && isvar('ri')
		yb = bi .* exp(-proj) + ri;
		yi = poisson(yb, [], 'try_poissrnd', 1); % kludge
		li_hat = -log((yi-ri) ./ bi);
		li_hat(yi-ri <= 0) = 0; % fix: need something better here...
	else
		li_hat = proj; % noiseless
	end
end


% FDK cone-beam reconstruction
if ~isvar('xfdk'), printm 'fdk'
	xfdk = feldkamp(cg, ig, li_hat, ...
		'extrapolate_t', ceil(1.3 * cg.nt/2)); % todo: compute carefully
%		'window', 'hanning,0.7', ... % test window
%	clf, im_toggle(xtrue, xfdk, [0 0.02]), return
%	im_toggle(permute(xtrue, [1 3 2]), permute(xfdk, [1 3 2]), [0 0.02])
prompt
end


if 0 % debugging tests
	if 1, printm 'nthread'
		tmp = feldkamp(cg, ig, li_hat, 'use_mex', 1, 'nthread', 1);
		jf_equal(xfdk, tmp)
	end

	if 1, printm 'mex2' % compare mex versions
		tmp = feldkamp(cg, ig, li_hat, 'use_mex', 2);
		max_percent_diff(xfdk, tmp)
	end

	if 1, printm 'mex3' % compare mex versions
		tmp = feldkamp(cg, ig, li_hat, 'use_mex', 3);
		max_percent_diff(xfdk, tmp)
	end

	if 1, printm 'no mex' % compare mex vs non-mex
		tmp = feldkamp(cg, ig, li_hat, 'use_mex', 0);
%		tmp = xfdk0 ~= 0; % circular mask
%		max_percent_diff(xfdk(tmp), xfdk0(tmp))
		max_percent_diff(xfdk, tmp)
%		clf, im_toggle(xfdk0, xfdk, [0 0.02]), return
	end
return
end


if im && ~isempty(xfdk)
	% show results (off-center slices worse than central slice)
	im(2, xfdk, 'FDK recon'), cbar
	im(3, xfdk - xtrue, 'FDK error'), cbar

	im subplot 5
	ix = 1:ig.nx; iy = ceil(ig.ny/2); iz = ceil(ig.nz/2);
	plot(ix, xtrue(ix,iy,iz), '-', ix, xfdk(ix,iy,iz), '--')
	axis([1 ig.nx -0.003 0.023])
%	legend('true', 'FDK recon', 'location', 'southoutside')
	title 'middle slice', xlabel 'ix'
	xtick(ix([1 end])), ytick([0 0.01 0.02])

	im subplot 6
	iz=1:ig.nz; ix = 1+floor(ig.nx/2); iy = 1+floor(ig.ny/2);
	plot(iz, squeeze(xtrue(ix,iy,iz)), '-', iz, squeeze(xfdk(ix,iy,iz)), '--')
	axis([1 ig.nz -0.003 0.023])
%	axis([1 ig.nz 0.007 0.013]) % cyl
%	legend('true', 'FDK recon', 2), xlabel 'iz'
	titlef('profile at (ix,iy)=(%g,%g)', ix,iy)
	xlabel iz
	xtick(iz([1 end])), ytick([0 0.01 0.02])
prompt
end

minmax(xtrue)
minmax(xfdk - xtrue, 'FDK error')
printm('nrmse = %g%%', nrms(xfdk(:), xtrue(:)) * 100)

if 0
	im plc
end
