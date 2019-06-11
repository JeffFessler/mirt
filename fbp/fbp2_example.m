% fbp2_example.m
% Example of how to use fbp2.m
% This compares "ordinary" FBP with FBP based on Mojette sampling.
% todo: this should eventually supercede fbp_fan_arc_example.m
% Copyright 2005-12-16, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'sino'
	down = 2;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500);
	ig = ig.downsample(down);

	% parallel-beam
	sg = sino_geom('par', 'nb', 888, 'na', 984, ...
		'orbit', 180, ... % usual case
		'strip_width', 'dr', ... % not important for FBP
		'dr', 541/949, 'offset_r', 0.25);
%		'orbit', 360, ... % 2008-10-14 just for testing
	sg = sg.downsample(down);

%	ell = [0 0 200 200 0 1; 0 0 10 10 0 0];
%	ell = [30 0 120 100 0 1; 0 0 10 10 0 0];
%	ell = [20 -15 225 150 30 10; 40 0 10 15 0 1];
	ell = []; clim = (1 + [-1 1] * 0.05) * 1000;
%	ell = [100 70 f.ds/2 f.ds/2 0 100]; % point source

	[xtrue ell] = ellipse_im(ig, ell, 'oversample', 4, 'rot', 90, ...
		'hu_scale', 1000);
	sino = ellipse_sino(sg, ell, 'oversample', 4);

	im plc 2 3
	im(1, xtrue, 'xtrue', clim), cbar
	im(4, sino, 'sino'), cbar

	t = floor(min(sg.nb * sg.d, ig.fov)/2);
	ig.mask = ellipse_im(ig, [0 0 t t 0 1]) > 0;
	im(6, ig.mask)
prompt
end


%
% conventional FBP reconstruction
%
if ~isvar('r_std'), printm 'fbp std'
	fg.std = fbp2(sg, ig);
	cpu etic
	r_std = fbp2(sino, fg.std);
	cpu etoc 'fbp std recon time'

	im(2, r_std, 'FBP matlab std', clim), cbar
	im(5, r_std - xtrue, 'error'), cbar
prompt
end

if 0 % compare single thread time
	fg.std1 = fbp2(sg, ig, 'nthread', 1);
	cpu etic
	r_std = fbp2(sino, fg.std1);
	cpu etoc 'fbp std recon time, nthread=1'
return
end


% dsc is too slow to be worth it
if 0 && ~isvar('r_dsc'), printm 'fbp dsc'
	fg.dsc = fbp2(sg, ig, 'type', 'dsc');
	cpu etic
	r_dsc = fbp2(sino, fg.dsc);
	cpu etoc 'fbp dsc recon time'

	im(3, r_dsc, 'FBP matlab dsc', clim), cbar
	im(6, r_dsc - xtrue, 'error'), cbar
prompt
end


%
% Mojette FBP reconstruction
% (an experimental approach - not essential)
%
if has_mex_jf
	if ~isvar('r_moj'), printm 'fbp moj'
		cpu etic
		fg.moj = fbp2(sg, ig, 'type', 'mojette', 'nthread', 1);
		cpu etoc 'moj setup time'

		cpu etic
		r_moj = fbp2(sino, fg.moj);
		cpu etoc 'moj recon time'

		im(3, r_moj, 'FBP mojette', clim), cbar
		im(6, r_moj - xtrue, 'error'), cbar
	prompt
	end

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
	if has_mex_jf
	im(3, r_moj, 'FBP mojette', clim), cbar
	im(6, r_moj - xtrue, 'error', elim), cbar
	end
end
