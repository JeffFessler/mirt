% fbp2_motion_example.m
% Example of effects of motion on fbp2.m
% Copyright 2007-10-12, Jeff Fessler, The University of Michigan

if ~isvar('xtrue'), printm 'xtrue'
	down = 4;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	sg = sino_geom('fan', 'nb', 888, 'na', 984, ...
		'dsd', 949, 'dod', 408, 'dfs', 0, ...
		'orbit', 720, ...
		'ds', 1, 'offset_s', 0.25, 'down', down);

	ell = []; clim = (1 + [-1 1] * 0.05) * 1000;

	[xtrue ell] = ellipse_im(ig, ell, 'rot', 90, 'hu_scale', 1000);
	ell0 = [ell; [-30 130 9 9 0 50]];
	xtrue0 = ellipse_im(ig, ell0, 'oversample', 4);

	ell1 = [ell; [-50 110 9 9 0 50]];
	xtrue1 = ellipse_im(ig, ell1, 'oversample', 4);

	im plc 2 3
	im(1, ig.x, ig.y, xtrue0, 'xtrue0', clim), cbar
	im(2, ig.x, ig.y, xtrue1, 'xtrue1', clim), cbar
	im(3, ig.x, ig.y, (xtrue0+xtrue1)/2, 'xtrue avg', clim), cbar

	t = floor(min(sg.nb * sg.d, ig.fov)/2);
	ig.mask = ellipse_im(ig, [0 0 t t 0 1]) > 0;
	im(6, ig.mask), cbar
prompt
end

%
% sinogram corrupted by motion
%
if ~isvar('sino'), printm 'sino'
	% moving ellipse strum
%	es = ellipse_motion(ell, 'type', 'none');
	es = ellipse_motion(ell0, 'type', 'linear', 'ellend', ell1);
	sino = ellipse_sino(sg, es, 'oversample', 4);
%	sinoo = ellipse_sino(sg, ell, 'oversample', 4);
%	minmax(sino2-sino)
	im(4, sino, 'sino'), cbar
end

%
% conventional FBP reconstruction
%
if ~isvar('r_std'), printm 'fbp std'
	fg.std = fbp2(sg, ig);
	cpu etic
	r_std = fbp2(sino, fg.std);
	cpu etoc 'fbp std recon time'

	im(5, r_std, 'FBP matlab std', clim), cbar
	im(6, r_std - xtrue, 'error'), cbar
prompt
end


return

max_percent_diff(xtrue, r_std)
max_percent_diff(xtrue, r_moj)
printm('sums: %g %g', [sum(r_std(:)) sum(r_moj(:))] / sum(xtrue(:)))
nrms(r_std, xtrue)
nrms(r_moj, xtrue)

if 0 % profiles
	ix=1:ig.nx; iy=1:ig.ny; ix=round(ig.nx*0.1973); ii=iy;
	clf, plot(ii, xtrue(ix,iy), 'y:', ...
		ii, r_std(ix,iy), 'c-', ii, r_moj(ix,iy), 'g--')
	axis([[0.14 0.7]*ig.nx clim])
	legend('true', 'FBP', 'Mojette', 'location', 'south')
return
end

if ~has_aspire, return, end % check consistency with aspire

if ~isvar('r_asp')
	dir = test_dir;
	f.sino = [dir 'sino.fld'];
	f.image = [dir 'image.fld'];
	f.dsc = [dir 't.dsc'];
	fld_write(f.sino, sino)

	t = aspire_pair(sg, ig);
	char_array_write(f.dsc, t)
	f.win = 'boxcar,1,0,1';
	com = sprintf('echo y | i -chat 0 fbp dsc %s %s %s %s', ...
		f.image, f.sino, f.dsc, f.win);
%	eval(['!' com])
	cpu etic
	disp(os_run(com))
	cpu etoc 'aspire time'

	if 0 % compare filters
		tmp = fld_read('fft_filt.fld');
		sum(tmp) / sum(test)
%		clf, plot([tmp test])
		max_percent_diff(test, tmp)
	return
	end

	if 0 % compare filtered projections
		% seem to match except for slight shift?? 
		tmp = fld_read('proj_filt.fld');
		clf, plot([tmp(:,1) test(:,1)])
		plot(tmp(:,1)-test(:,1))
		max_percent_diff(test, tmp)
		sum(tmp(:)) / sum(test(:))
		minmax(test)
		minmax(tmp)
		minmax(test-tmp)
	return
	end

	r_asp = fld_read(f.image);
	im(1, r_asp, 'aspire', clim), cbar
	max_percent_diff(r_std, r_asp)
prompt
end

if 1 % show all images
	im plc 2 3; clim = [0.99 1.04]*1000; elim = [-0.1 0.1]*1000;
	im(1, r_asp, 'FBP aspire', clim), cbar
	im(4, r_asp - xtrue, 'error', elim), cbar
	im(2, r_std, 'FBP matlab', clim), cbar
	im(5, r_std - xtrue, 'error', elim), cbar
	im(3, r_moj, 'FBP mojette', clim), cbar
	im(6, r_moj - xtrue, 'error', elim), cbar
end
