% Gtomo2_sparse_test.m
% test the Gtomo2_sparse object

if ~isvar('G'),	disp 'setup Gtomo2_sparse_test'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');

	if ~exist(f.dir, 'dir'), error 'need to create directory', end
	os_run(sprintf('wt -chat 0 dsc 2 nx 16 ny 14 nb 20 na 10 >! %s', f.dsc'))
	os_run(sprintf('echo y | wt -chat 0 gen %s', f.dsc'))

	if 1
		Gd = Gtomo2_dsc(f.dsc);
		na = Gd.na;
		nb = Gd.nb;
		Gm = Gd(:,:);
		G = Gtomo2_sparse(Gm, 1, Gd.nx, Gd.ny, nb, na, Gd.mask);
	else
		G = Gtomo2_sparse(f.wtf);
	end
	im(331, G.mask, 'mask')
prompt
end

if 1, disp 'test block object'
	nblock = 4;
	Gb = Gblock(G, nblock, 1);
	x = G.mask;
	istart = 3;
	ia = istart:nblock:na;
	yb = zeros(nb, na);
	yb(:,ia) = reshape(Gb{3} * x(:), nb, length(ia));
	ya = reshape(G * x(:), nb, na);

	if 1
		G = make_block(G, nblock);
		mtimes_block(G, x(:), istart, nblock);
	end
prompt
end

if 1
	x = G.mask;
	y = reshape(G * x(:), G.nb, G.na);
	im(332, y)

	y = ones(G.nb, G.na);
	x = reshape(G' * y(:), G.nx, G.ny);
	im(333, x)

	j = sub2ind([G.nx G.ny], G.nx/2, G.ny/2+2);
	t = reshape(G(:,j), G.nb, G.na);
	im(334, t)

	t = reshape(G(:,[j j+1]), [G.nb G.na 2]);
	im(335, t)

	i = sub2ind([G.nb G.na], G.nb/2, G.na/2);
	t = reshape(G([i i+2],:), [G.nx G.ny 2]);
	im(336, t)
prompt
end

if 1,	disp 'test adjoint'
	test_adjoint(G);
end
