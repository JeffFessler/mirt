% ir_fbp2_demo1.m
% demonstrate how FBP works with view-by-view backprojection
% Copyright 2016-03-07, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'sino'
	down = 2;
	ig = image_geom('nx', 504, 'ny', 512, 'fov', 500);
	ig = ig.downsample(down);

	% parallel-beam
	sg = sino_geom('par', 'nb', 888, 'na', 984, ...
		'orbit', 180, ... % usual case
		'strip_width', 'dr', ... % not important for FBP
		'dr', 541/949, 'offset_r', 0.25);
	sg = sg.downsample(down);

	[xtrue ell] = ellipse_im(ig, [], 'oversample', 2, ...
		'hu_scale', 1000);
	sino = ellipse_sino(sg, ell, 'oversample', 2);

	t = floor(min(sg.nb * sg.d, ig.fov)/2);
	ig.mask = ellipse_im(ig, [0 0 t t 0 1]) > 0;

	fstruct = fbp2(sg, ig);
end

if 1
	im plc 3 4
	clim = [900 1100];
	im(1, ig.x, ig.y, xtrue, 'Image', clim), cbar
	xtick([-1 0 1] * 240), ytick([-1 0 1] * 250)
	xlabel 'x [mm]'
	ylabel 'y'

	im(5, sg.s, sg.ad, sino, 'Sinogram'), cbar
	xlabel 'r [mm]'
	ylabel 'angle [degrees]'
	xtick([-1 0 1] * 250)
	axisy([0 180]), ytick([0 180])
end

view_set_list = {55, 1:12:sg.na, 1:sg.na};

[~, sino_filt] = fbp2(sino, fstruct);

for ivs = 1:numel(view_set_list)
	ia = view_set_list{ivs};
	tmp = sg.zeros;
	tmp(:,ia) = sino(:,ia) * sg.na / numel(ia); % scale to # of views
	out = fbp2(tmp, fstruct);
%	out = fbp2_back(sg, ig, tmp); % un-filtered

	if numel(ia) >= sg.na / 12
		tlim = clim;
	else
	%	tlim = minmax(out)';
		tlim = mean(col(out(end/2+[-11:11],end/2+[-11:11])));
		tlim = tlim + [-1 1] * 300;
	end
	im(ivs+1, ig.x, ig.y, out, tlim), cbar
	xtick([-1 0 1] * 240), ytick([-1 0 1] * 250)

	titlef('FBP of %d/%d views', numel(ia), sg.na)
	im(5 + ivs, sg.s, sg.ad, sino_filt, 'Filtered Sinogram'), cbar
	xtick([-1 0 1] * 250)
	axisy([0 180]), ytick([0 180])
	if numel(ia) < sg.na
		hold on
		r = sg.s;
		plot([r(1) r(sg.nb)], sg.ad(ia) * [1 1], 'b-')
		hold off
	end

	drawnow
end


if 0 % movie of FBP reconstruction
	for ia=1:sg.na
		tmp = sino;
		tmp(:,ia+1:end) = 0;
		out = fbp2(tmp, fstruct);

		tmp = sprintf('%d views', ia);
		im(3, out, tmp), cbar
		drawnow
	end
end
