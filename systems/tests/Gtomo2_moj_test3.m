% Gtomo2_moj_test3
% Test the Mojette parallel-beam backprojector vs std pixel-driven.
% Copyright 2006-1-31, Jeff Fessler, The University of Michigan

redo = 0;
im pl 2 3
if ~isvar('G'), printm 'G'
	down = 4;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 512, 'down', down);
%	ig.mask = ellipse_im(ig, [0 0 250*[1 1] 0 1], 'oversample', 2) > 0;

	f.table = {'mojette,back1'};

%	f.nr = ceil(sqrt(2)/2 * max(ig.nx,ig.ny))*2; % circle mask
%	f.nr = 2 * max(ig.nx,ig.ny) + 4; % square mask
	f.nr = 64+ceil(ig.fov / ig.dx * 1.5 / 2) * 2; % bug?
	sg = sino_geom('moj', 'nb', f.nr, 'na', 984/down, ...
		'dx', ig.dx, 'offset_r', 0, ...
		'orbit', 180, 'orbit_start', 0);

	cpu etic
	G = Gtomo2_table(sg, ig, f.table, 'nthread', 1);
	cpu etoc 'Gmoj precompute time:'
prompt
end

if 0
	sino = G * ig.ones;
	im(sino), cbar
return
end

if 1
%	sino = sg.ones;
%	sino = ones(sg.nb,1) * unitv(sg.na, 7)';
%	tt = zeros(sg.nb, 2);
	for ib=84%1:10:sg.nb
		ticker(mfilename, ib, sg.nb)
		sino = sg.unitv(ib, sg.na/2+1);
		x1 = G' * sino; % backproject

		arg1 = {uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, ig.offset_y};
		arg2 = {sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start, ...
			int32(1)}; % nthread
		x2 = jf_mex('back2', arg1{:}, arg2{:}, single(sino));
		x2 = x2 / single(pi) * sg.na;
		x2 = x2 .* ig.mask;

		im(1, x1), im(2, x2), im(3, x2-x1)
%		tt(ib,:) = minmax(x2-x1);
	end
%clf, plot(tt)
return

	for ia=1:sg.na
		ticker(mfilename, ia, sg.na)

		sino = ones(sg.nb,1) * unitv(sg.na, ia)';
		sino = sg.unitv(sg.nb/2, ia);
		x1 = G' * sino; % backproject

		arg2 = {sg.d(ia), sg.offset, sg.orbit, sg.ad(ia)};
		x2 = 0*x2 + jf_mex('back2', arg1{:}, arg2{:}, single(sino(:,ia)));
		x2 = x2 / single(pi); % * sg.d(ia) / ig.dx;
tt(ia,:) = minmax(x2-x1);
%max_percent_diff(x2,x1)
%		im(1, x1), im(2, x2), im(3, x2-x1)
%pause
	end
clf, plot(tt)
return
	x2 = x2 / single(pi);
	x2 = x2 .* ig.mask;

	im pl 1 3
	im(1, x1), im(2, x2), im(3, x2-x1)
end
