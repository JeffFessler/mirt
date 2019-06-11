% Cdiff_test.m
% test Cdiff object

if 1 % look at 2nd order case
	nx = 8; ny = 6;
	C = Cdiff(ones(nx,ny), 'order', 2);
	Cf = C(:,:);
	t = reshape(Cf', nx*ny, nx*ny, []);
	im(t)
	t = zeros(nx,ny); t(nx/2+1,ny/2+1) = 1;
	t = C' * (C * t);
	im(t)
end

if 1
	nx = 16; ny = 14;
	nx = 512; ny = 500; % for large images, the mex file is much faster!
	ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
	ig.mask = [0 0 [nx ny]/2-5 0 1];
	ig.mask = conv2(double(ellipse_im(ig, []) > 0), ones(2), 'same') > 0;
%	mask = ones(nx,ny); % all

	if 1, printm 'test penalty_mex'
		x = single(ig.mask);
		offsets = [nx-1];
		offsets = int32(offsets);
		d1 = penalty_mex('diff2,forw1', x, offsets);
		d2 = penalty_mex('diff2,forw1', x, offsets, int32(ndims(x)));
		if any(d1(:) ~= d2(:)), error 'bug', end

		x1 = penalty_mex('diff2,back1', d1, offsets);
		x2 = penalty_mex('diff2,back1', d1, offsets, int32(ndims(x)));
		if any(x1(:) ~= x2(:)), error 'bug', end
	end

	ctype = 'leak';
	ctype = 'tight';
	order = 1;
	tic
	C1 = Cdiff(ig.mask, 'edge_type', ctype, 'offsets', '2d,hvd', ...
		'distance_power', 1., 'order', 1);
	printm('make C1 time %g', toc)

	tic
	[C2 wjk] = C2sparse(ctype, ig.mask, 8);
	printm('make C2 time %g', toc)
	C2 = spdiag(sqrt(wjk)) * C2;
	C2 = C2(:,ig.mask(:));

	if 0 % test old Cmask
		cpu etic
		b1 = Cmask('tight,2d,hvd', ig.mask);
		cpu etoc 'make scale time:'
		b1(:,:,[3 4]) = b1(:,:,[3 4]) / sqrt(sqrt(2));
		if 0
			b1 = reshape(sqrt(wjk), [nx ny 4]);
		end

		cpu etic
		b2 = penalty_mex('scales,tight', single(ig.mask), C1.arg.offsets, 1.);
		b2 = double(b2);
		cpu doc 'make scale time:'

		im plc 1 3, im(1, b1), im(2, b2), im(3, b1-b2)
		printm('old vs new: %g%%', max_percent_diff(b1, b2))
		equivs(b1, b2)
	return
	end

	rng(0)
	x = rand(nx, ny);
	x = dsingle(x);
	x = x .* ig.mask;

%prompt
end

if 1
	xm = double(x(ig.mask(:)));
	cpu etic
	d1 = C1 * x; 
	cpu etoc 'C1 forw time:'

	cpu etic
	d2 = C2 * xm; 
	cpu etoc 'C2 forw time:'
	d2 = reshape(d2, [nx ny 4]); 
%	printm('Cx vs penalty_mex: %g%%', max_percent_diff(d1, d2))
	equivs(d1, d2)

	if im
		im plc 1 3
		im(1, d1, 'C * x'), cbar
		im(2, d2, 'penalty--mex'), cbar
		im(3, d2-d1, 'err'), cbar
	prompt
	end
end

if 1 % [x x]
	xx = double([xm xm]);
	cpu etic
	d11 = C1 * xx;
	cpu etoc 'C1 forw time:'

	cpu etic
	d22 = C2 * xx; 
	cpu etoc 'C2 forw time:'
%	printm('Cxx vs penalty_mex: %g%%', max_percent_diff(d11, d22))
	equivs(d11, d22)
end

if 1
	d = double(d1(:));
	%d = zeros(nx,ny,4);
	%d(end/2,end/2,1) = 1;

	cpu etic
	x1 = C1' * d;
	cpu etoc 'C1''d time:'
	x1 = embed(x1, ig.mask);

	cpu etic
	x2 = C2' * d;
	cpu etoc 'C2''d time:'
	x2 = embed(x2, ig.mask);

%	printm('C1 vs C2: %g%%', max_percent_diff(x1, x2))
	equivs(x1, x2)

	if im
		im plc 1 3
		im(1, x1, 'C''d'), cbar
		im(2, x2, 'penalty--mex'), cbar
		im(3, x2-x1, 'err'), cbar
	end
end

if 1 % C' [d d]
	dd = double([d(:) d(:)]);

	cpu etic
	x11 = C1' * dd;
	cpu etoc 'C1''d time:'

	cpu etic
	x22 = C2' * dd;
	cpu etoc 'C2''d time:'

%	printm('C1 vs C2 for dd: %g%%', max_percent_diff(x11, x22))
	equivs(x11, x22)
end
