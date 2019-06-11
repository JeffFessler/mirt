% fbp_fan_arc_example.m
% example of how to use fbp_fan_arc.m
% Copyright Jeff Fessler, University of Michigan

if ~isvar('sino')
	down = 4;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	sg = sino_geom('fan', 'nb', 888, 'na', 984, 'ds', 1.0, ...
		'dsd', 949, 'dod', 408, 'offset_s', 0.25, ...
		'orbit', 360, ...
		'down', down);

	ell = [0 0 200 200 0 1; 0 0 10 10 0 1];
	ell = [0 0 200 200 0 0; 0 0 10 10 0 1];
	ell = [20 -15 225 150 30 10; 40 0 10 15 0 1];
%	ell = [100 100 ds/4 ds/4 0 1]; % point source

	sino = ellipse_sino(sg, ell, 'oversample', 2);

	xtrue = ellipse_im(ig, ell, 'oversample', 4);

	% system object
	G = Gtomo2_dscmex(sg, ig);

	im plc 2 3
	im(1, xtrue, 'x'), cbar
	im(4, sino, 'sino'), cbar
prompt
end

if 0 % examine recon of flat fan-beam sinogram
	sino = repmat(abs(sg.s) < 0.8 * max(sg.s), [1 sg.na]);
	im(sino)
end

% fan-beam reconstruction
if ~isvar('recon')
	recon = fbp_fan_arc(sino, G, 'ramp');
end

if 1 % compare to new fbp2() method
	tmat = fbp2(sg, ig, 'type', 'std:mat');
	rmat = fbp2(sino, tmat);
	max_percent_diff(rmat, recon)
prompt
end

if 0 % compare to new fbp2() method
	tmex = fbp2(sg, ig, 'type', 'std:mex');
	rmex = fbp2(sino, tmex);
	mmask = (rmat > 0);
	max_percent_diff(rmex .* mmask, recon .* mmask)
	im(2, rmat .* mmask, 'mat'), cbar
	im(3, rmex .* mmask, 'mex'), cbar
return
end

if 0 % examine "aliasing"
	clim = [9.5 10.5];
%	clim = [0.8 1.2];
%	clim = [-0.2 0.2];
	im clf, im(recon, clim), cbar
return
end

if im
	im(2, recon, 'FBP matlab'), cbar
	im(5, recon - xtrue, 'error'), cbar
	iy = ig.ny/2; ix = 1:ig.nx;
	subplot(133)
	plot(ix, xtrue(ix,iy), '-', ix, recon(ix,iy), '--')
	axis([1 ig.nx -0.5 10.5]), legend('true', 'recon', 4)
end

if has_aspire % check consistency with aspire
	dir = test_dir;
	f.sino = [dir 'sino.fld'];
	f.image = [dir 'image.fld'];
	f.dsc = [dir 't.dsc'];
	fld_write(f.sino, sino)

	char_array_write(f.dsc, G.arg.args)
	f.win = 'boxcar,1,0,1';
	com = sprintf('echo y | i fbp dsc %s %s %s %s', ...
		f.image, f.sino, f.dsc, f.win)
%	eval(['!' com])
	disp(os_run(com))

	if 0 % compare filters
		tmp = fld_read('fft_filt.fld');
		sum(tmp) / sum(test)
%		im clf, plot([tmp test])
		max_percent_diff(test, tmp)
	return
	end

	if 0 % compare filtered projections
		% seem to match except for slight shift?? 
		tmp = fld_read('proj_filt.fld');
		im clf, plot([tmp(:,1) test(:,1)])
		plot(tmp(:,1)-test(:,1))
		max_percent_diff(test, tmp)
		sum(tmp(:)) / sum(test(:))
		minmax(test)
		minmax(tmp)
		minmax(test-tmp)
	return
	end

	im_asp = fld_read(f.image);
	im(233, im_asp, 'aspire'), cbar
	good = recon ~= 0;
	diff = im_asp - recon;
	diff = diff .* good;
	im(236, diff, 'aspire-matlab'), cbar
	im(234, (im_asp ~= 0) - (recon ~= 0), 'aspire-matlab support'), cbar
	max_percent_diff(recon.*good, im_asp.*good)
end
