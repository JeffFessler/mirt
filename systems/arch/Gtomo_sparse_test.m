% Gtomo_sparse_test.m
% test the Gtomo_sparse object

if ~has_aspire
	disp 'skipping since no aspire'
	return
end

if ~isvar('Gm'), disp 'setup Gtomo_sparse_test'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');

	if ~exist(f.dir, 'dir'), error 'need to create directory', end
	os_run(sprintf('wt -chat 0 dsc 2 nx 16 ny 14 nb 20 na 10 >! %s', f.dsc'))
	os_run(sprintf('echo y | wt -chat 0 gen %s', f.dsc'))

	Gd = Gtomo2_dsc(f.dsc);
	Gm = Gd(:,:);
prompt
end

if ~isvar('G'),	disp 'G'
	if 1
		nx = Gd.nx;
		ny = Gd.ny;
		na = Gd.na;
		nb = Gd.nb;
		G = Gtomo_sparse(Gm, Gd.mask, 'idim', [nx ny], 'odim', [nb na]);
	else
		G = Gtomo_sparse(f.wtf);
		nx = G.arg.idim(1);
		ny = G.arg.idim(2);
		nb = G.arg.odim(1);
		na = G.arg.odim(2);
	end
	im clf
	im(331, G.arg.mask, 'mask')
prompt
end

if 1
	x = G.arg.mask;
	y = G * x;
	im(332, y, 'G*x'), cbar

	y = ones(nb, na);
	x = G' * y;
	im(333, x, 'G''*1'), cbar

	j = find(G.arg.mask);
	j = find(j == sub2ind([nx ny], nx/2, ny/2+2))
	t = reshape(G(:,j), nb, na);
	im(334, t, 'G(:,j)'), cbar

	t = reshape(G(:,[j j+1]), [nb na 2]);
	im(335, t, 'G(:,[j j+1])')

	i = sub2ind([nb na], nb/2, na/2);
	t = embed(G([i i+2],:), G.arg.mask);
	im(336, t, 'G([i i+1],:)')
prompt
end

if 1,	disp 'test adjoint'
	test_adjoint(G);
end

if 1, disp 'test block object'
disp 'todo: WARNING - BLOCK TESTING NOT DONE!'
return
	nblock = 4;
	Gb = Gblock(G, nblock, 1);
	x = G.arg.mask;
	istart = 3;
	ia = istart:nblock:na;
	yb = zeros(nb, na);
	yb(:,ia) = reshape(Gb{3} * x(:), nb, length(ia));
	ya = reshape(G * x(:), nb, na);
prompt
end
