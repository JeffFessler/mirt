% Gblock_test.m
% Test the Gblock object in all of its versions.
% This also tests several system objects:
% Gsparse, Gtomo2_dscmex, and Gtomo2_wtmex
%
% Copyright 2002-2-19, Jeff Fessler, University of Michigan

% image and sparse matrix
if ~isvar('A0'), printm 'setup Gblock_test'
	ig = image_geom('nx', 20, 'ny', 22, 'dx', 1);
	sg = sino_geom('par', 'nb', 24, 'na', 18, ...
		'ray_spacing', 1, 'strip_width', 1);
	xt = ellipse_im(ig, 'shepplogan-emis');
	ig.mask = ig.circ > 0;

	f.strip_width = 1;
	f.strips = {'strip_width', f.strip_width};
	A1 = Gtomo2_strip(sg, ig, f.strips{:});
%	A0 = A1.arg.G; % [nd np]
	A0 = A1.arg.matrix; % [nd np]

	im plc 3 3
	im(1, ig.mask, 'mask')
	im(2, xt, 'xt')

	xm = xt(ig.mask);
	ym = sg.shape(A0 * double(xm)); % matlab sparse needs double
	yt = A1 * xt;
	im(3, ym)
	jf_equal(yt, A1 * xt)
prompt
end


% 2d system objects
if ~isvar('As'), printm 'As'
	As = Gsparse(A0, 'mask', ig.mask, 'odim', [sg.nb sg.na]);
end

if has_aspire && ~isvar('Aw'), printm 'Gtomo2_wtmex'
	Aw = Gtomo2_wtmex(sg, ig, 'grouped', 'row', ...
		'pairs', {'tiny', 0, f.strips{:}});
	yy = Aw * xt;
	equivs(yy, yt)
	clear yy
prompt
end

if has_aspire && ~isvar('Ad'), printm 'Gtomo2_dscmex'
	f.dir = test_dir;
	f.mask = [f.dir 'mask.fld'];
	fld_write(f.mask, ig.mask, 'check', 0)

	if 0
		f.arg = arg_pair('system', 2, ...
			'nx', nx, 'ny', ny, 'nb', nb, 'na', na, ...
			'orbit', 180, 'orbit_start', 0, ...
			'pixel_size', 1, 'ray_spacing', 1, 'strip_width', 1, ...
			'scale', 0);
%		Ad = Gtomo2_dscmex(f.arg);
	end

	Ad = Gtomo2_dscmex(sg, ig);

	yy = Ad * xt;
	equivs(yy, yt)
	clear yy
prompt
end


% basic tests/comparisons of the objects
if 1
	if has_aspire
		Alist = {A1, As, Aw, Ad};
	else
		Alist = {A1, As};
	end

	for ii = 1:length(Alist)
		A = Alist{ii};
		pr class(A)

		% check sum
		t1 = sum(A1);
		t2 = sum(A);
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('sum	mpd %g', mpd)
			if mpd/100 > 1e-6, error sum, end
		end

		% check A*x
		t1 = A1 * xm;
		t2 = A * xm;
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A*x	mpd %g', mpd)
			if mpd/100 > 1e-6, error Ax, end
		end

		% check A*[x x]
		t1 = A1 * [xm xm];
		t2 = A * [xm xm];
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A*[x x]	mpd %g', mpd)
			if mpd/100 > 1e-6, error Axx, end
		end

		% check A'*y
		t1 = A1' * yt;
		t2 = A' * yt;
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A''y	mpd %g', mpd)
			if mpd/100 > 1e-6, error Aty, end
		end

		% check A'*[y y]
		t1 = A1' * [yt(:) yt(:)];
		t2 = A' * [yt(:) yt(:)];
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A''[y y]	mpd %g', mpd)
			if mpd/100 > 1e-6, error Atyy, end
		end

		% check A(:,j)
		j = find(ig.unitv);
		t1 = full(A1(:,j));
		t2 = full(A(:,j));
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A(:,j)	mpd %g', mpd)
			if mpd/100 > 1e-6, error A(:,j), end
		end

		% check A(:,js)
		j = [0 1] + j;
		t1 = full(A1(:,j));
		t2 = full(A(:,j));
		equivs(t1, t2)
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('A(:,js)	mpd %g', mpd)
			if mpd/100 > 1e-6, error A(:,js), end
		end

		% check A(:,:)
		if ig.nx < 100
			t1 = A1(:,:);
			t2 = A(:,:);
			equivs(t1, t2)
			if 0
				mpd = max_percent_diff(t1,t2);
				printf('A(:,:)	mpd %g', mpd)
				if mpd/100 > 1e-6, error A(:,:), end
			end
		end
	end

prompt
end


% block versions of each object
if ~isvar('Bs'), printm 'setup Gblock_test'
	nblock = 8;
	Bs = Gblock(As, nblock);
	if has_aspire
		Bd = Gblock(Ad, nblock);
		Bw = Gblock(Aw, nblock);
	end
prompt
end


% check acceleration
if 0
%	profile on
	cpu etic
	for ii=1:10, y1 = As * xm; end
	cpu etoc 'As time:'
	cpu etic
	for ii=1:10, y2 = Bs{1} * xm; end
	cpu etoc 'Ab time:'
%	profile report
return
end


% basic block tests
if 1, printm 'block tests'
	if has_aspire
		Blist = {Bs, Bd, Bw};
	else
		Blist = {Bs};
	end

	for ii = 1:length(Blist)
		B = Blist{ii};

		% check mask
		if 1
			t1 = ig.mask;
			t2 = B.arg.mask;
			if isempty(t2), t2 = B.arg.ig.mask; end
			jf_equal(t1, t2)
		end

		% check Bx
		t1 = A1 * xm;
		t2 = B * xm;
		equivs(t1,t2);
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('B*x	mpd %g', mpd)
			if mpd/100 > 1e-6, error Bx, end
		end

		% check B'y
		t1 = A1' * yt;
		t2 = B' * yt;
		equivs(t1,t2);
		if 0
			mpd = max_percent_diff(t1,t2);
			printf('B''*y	mpd %g', mpd)
			if mpd/100 > 1e-6, error Bty, end
		end

		% block operations
		for k=1:nblock
			ia = k:nblock:sg.na;

			% check B{k}*x
			t1 = A1 * xm;
			t1 = sg.shape(t1);
			t1 = t1(:,ia);
			t1 = t1(:);
			t2 = B{k} * xm;
			equivs(t1,t2)
			if 0
				mpd = max_percent_diff(t1,t2);
				printf('B{k}*x	mpd %g', mpd)
				if mpd/100 > 1e-6, error B{k}*x, end
			end

			% check B{k}'*y()
			tmp = sg.zeros;
			tmp(:,ia) = yt(:,ia);
			t1 = A1' * tmp(:);
			t2 = B{k}' * col(yt(:,ia));
			equivs(t1,t2)
			if 0
				mpd = max_percent_diff(t1,t2);
				printf('B{k}''*y	mpd %g', mpd)
				if mpd/100 > 1e-6, error 'B{k}''*y()', end
			end

		end
	end
prompt
end
