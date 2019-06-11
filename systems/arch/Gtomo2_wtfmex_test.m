% Gtomo2_wtfmex_test.m
% Test the Gtomo2_wtfmex object

if ~has_aspire
	return
end

% test new vewsion
if 1
	ig = image_geom('nx', 32, 'ny', 30, 'dx', 2);
	sg = sino_geom('par', 'nb', 40, 'na', 20, 'dr', 1.8);
	ig.mask = ig.circ(ig.fov/2) > 0;
	im(ig.mask)
	G = Gtomo2_wtmex(sg, ig, 'chat', 0, 'nthread', 2);
	Gtomo2_test(G, ig.mask) % put it through paces
	clear G
end

if ~isvar('G'),	printm 'setup Gtomo2_wtfmex_test'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');

	os_run(sprintf('wt -chat 0 dsc 2 >! %s', f.dsc))
	os_run(sprintf('echo y | wt -chat 0 gen %s', f.dsc))

	G = Gtomo2_wtfmex(f.wtf);
	im(331, G.mask, 'mask')
prompt
end

if 1,	printm 'basic tests'
	x = G.mask;
	y = reshape(G * x(:), G.nb, G.na);
	im(332, y)

	y = ones(G.nb, G.na);
	x = reshape(G' * y(:), G.nx, G.ny);
	im(333, x, 'x'), cbar

	x2 = reshape((G.^2)' * y(:), G.nx, G.ny);
	im(337, x2, 'x2'), cbar

	t = reshape(G(:,2718), G.nb, G.na);
	im(334, t)

	t = reshape(G(:,[2718 2710]), [G.nb G.na 2]);
	im(335, t)

	t = reshape(G([202 210],:), [G.nx G.ny 2]);
	im(336, t)
prompt
end

if 1, printm 'adjoint test'

	os_run(sprintf('wt -chat 0 dsc 2 nx 12 ny 14 nb 16 na 10 >! %s', f.dsc))
	os_run(sprintf('echo y | wt -chat 0 gen %s', f.dsc))
	G = Gtomo2_wtfmex(f.wtf);
	test_adjoint(G);
	clear G
end
