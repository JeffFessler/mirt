% fbp_fan_arc_point
% test FBP of point-like object with various offsets

if ~isvar('sino')
	down = 2;
	ig = image_geom('nx', 128, 'ny', 120, 'fov', 30, 'down', down);
	sg = sino_geom('fan', 'nb', 888, 'na', 984, 'ds', 1.0, ...
		'dsd', 949, 'dod', 408, 'offset_s', 0.25, ...
		'down', down);

	ell = [2 1 0.5 0.5 0 1];
	sino = ellipse_sino(sg, ell, 'oversample', 3);

	x = ellipse_im(ig, ell, 'oversample', 4);

	G = Gtomo2_dscmex(sg, ig); % system object

	im plc 2 3
	im(1, x, 'x'), cbar
	im(2, sino, 'sino'), cbar
prompt
end


if ~isvar('xs'), printm 'try various offsets'
	offsets = [-0.75:0.25:1];

	xs = ig.zeros('nz', length(offsets));

	for io=1:length(offsets)

		offset = offsets(io);
		ss = sg;
		ss.offset_s = offset;
		G = Gtomo2_dscmex(ss, ig);

		xs(:,:,io) = fbp_fan_arc(sino, G);
		im clf, im(xs(:,:,io))
		if im, title(sprintf('offset = %g', offset)), drawnow, end
	end
prompt
end
im(xs)

im plc 3 4
no = length(offsets);
for io=1:no
	x = xs(:,:,io);
	im(io, x, sprintf('offset = %g', offsets(io)))
	fw(io) = fwhm2(x-min(x(:)), 'dx', ig.dx);
end
if im
	subplot(313), plot(offsets, fw, '-o')
	xtick(offsets), axis([minmax(offsets)' 0.5 2.5]), grid
	xlabel 'channel offset', ylabel 'FWHM [mm]'
	title 'fbp-fan-arc-point'
end
