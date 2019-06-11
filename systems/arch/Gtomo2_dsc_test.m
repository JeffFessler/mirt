% Gtomo2_dsc_test.m
% test the Gtomo2_dsc object

if ~has_aspire
	return
end

if ~isvar('G'),	disp 'setup Gtomo2_dsc_test'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.wtf = '';		% use this to avoid testing wtf version
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');
	f.nx = 64;
	f.ny = 62;
	f.na = 30;
	f.system = 2;

	f.com = sprintf('wt -chat 0 dsc %d nx %d ny %d na %d scale 0 >! %s', ...
		f.system, f.nx, f.ny, f.na, f.dsc);
	os_run(f.com)

	if f.wtf
%		if exist(f.wtf), delete(f.wtf), end
		os_run(sprintf('echo y | wt -chat 0 gen %s', f.dsc))
	end

	G = Gtomo2_dsc(f.dsc, 1, 0);
	na = G.na;
	nb = G.nb;
	nx = G.nx;
	ny = G.ny;
	im clf, im(331, G.mask, 'Gtomo2\_dsc mask')
prompt
end

if ~isvar('A') & f.wtf
	A = wtfmex('load', f.wtf);
end

%
% timing with threads
%
if 1,	disp 'thread timing test'
	if ~isvar('G2')
		G2 = Gtomo2_dsc(f.dsc, 2, 0);
	end

	tic
	disp 'starting threaded one'
	for ii=1:4, t = G2 * ones(nx*ny,1); end
	printf('threaded time = %g', toc)

	tic
	disp('starting ordinary one')
	for ii=1:4, t = G * ones(nx*ny,1); end
	printf('non-threaded time = %g', toc)
prompt
end

%
% test basic stuff
%
if 1
	t = reshape(G(:,[2718 2710]), [nb na 2]);
	im(332, t, 'columns')
	t = reshape(G([202 310],:), [nx ny 2]);
	im(333, t, 'rows')

	x = double(G.mask);
	y1 = reshape(G * x(:), nb, na);
	if isvar('A')
		y2 = reshape(A * x(:), nb, na);
	else
		y2 = zeros(size(y1));
	end
	im(334, y1, 'forward G'), cbar
	im(335, y2, 'forward A'), cbar
	im(336, y2-y1, 'A-G'), cbar
	printf('forward error = %g%%', max_percent_diff(y1,y2))

	y = ones(nb, na);
	x1 = reshape(G' * y(:), nx, ny);
	if isvar('A')
		x2 = reshape(A' * y(:), nx, ny);
	else
		x2 = zeros(size(x1));
	end
	im(337, x1, 'backward G'), cbar
	im(338, x2, 'backward A'), cbar
	im(339, x2-x1, 'A-G'), cbar
	printf('back error = %g%%', max_percent_diff(x1,x2))
end

%
% test angular subsets aspects
%
if 0
	ia = 1:6:na;
	ii = outer_sum([1:nb]',(ia-1)*nb);

	y = ones(nb,na);
	Gt = G';
	x = reshape(Gt(:,ii) * y(ii), nx, ny);
	im clf, im(221, x, 'subset'), cbar

	x = double(G.mask);
	y1 = zeros(nb,na);
	t = Gt(:,ii)';
	y1(:,ia) = reshape(Gt(:,ii)' * x, nb, length(ia));
	y2 = zeros(nb,na);
	y2(:,ia) = reshape(A(ii,:) * x(:), nb, length(ia));
	im(222, y1, 'subset G'), cbar
	im(223, y2, 'subset A'), cbar
	im(224, y2-y1, 'A-G'), cbar
prompt
end

%
% test adjoint
%

if ~isvar('Gs'),	disp 'setup for adjoint'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.nx = 16;
	f.ny = 14;
	f.na = 12;

	f.com = sprintf('wt -chat 0 dsc %d nx %d ny %d na %d scale 0 >! %s', ...
		f.system, f.nx, f.ny, f.na, f.dsc);
	os_run(f.com)
	Gs = Gtomo2_dsc(f.dsc, 1, 0);
prompt
end

if 1, disp 'adjoint test'
	test_adjoint(Gs);
end
