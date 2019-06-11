% time_Gsparse.m
% timing test for Gsparse vs Gtomo2_sparse

if ~isvar('Gs'), disp 'Gs'
	f.dir = test_dir;
	f.dsc = [f.dir 't.dsc'];
	f.wtf = strrep(f.dsc, 'dsc', 'wtf');
	os_run(['wt -chat 0 dsc 12 nx 64 fwhm_detector 5 >! ' f.dsc])
	os_run(['echo y | wt -chat 0 gen ' f.dsc])
	Gs = Gsparse(f.wtf);
end

if ~isvar('G2b'), disp 'G2b'
	f.nblock = 6;
	Gsb = Gblock(Gs, f.nblock, 0);
	nx = Gs.arg.idim(1);
	ny = Gs.arg.idim(2);
	nb = Gs.arg.odim(1);
	na = Gs.arg.odim(2);
	mask = Gs.arg.mask;
	G2 = Gtomo2_sparse(Gs.arg.G, 0, nx, ny, nb, na, mask);
	G2b = Gblock(G2, f.nblock, 0);
end

x = double(mask(mask));
xs = single(mask(mask));
if 1, disp 'time'
	y = Gs * xs;
	y = G2 * x;
	x = Gs' * y;
	x = G2' * y; % warm-up
	cpu etic, for ii=1:3, y = Gs * x; end, cpu etoc 'Gs time'
	cpu etic, for ii=1:3, y = G2 * x; end, cpu etoc 'G2 time'
	ys = single(y);
	cpu etic, for ii=1:3, x = Gs' * y; end, cpu etoc 'Gs time'
	cpu etic, for ii=1:3, x = G2' * y; end, cpu etoc 'G2 time'
end

if 1
	y = Gsb{1} * x;
	y = G2b{1} * x;
	x = Gsb{1}' * y;
	x = G2b{1}' * y; % warm-up

	cpu etic, for ii=1:f.nblock, y = Gsb{ii} * x; end, cpu etoc 'Gs time'
	cpu etic, for ii=1:f.nblock, y = G2b{ii} * x; end, cpu etoc 'G2 time'
	ys = single(y);
	cpu etic, for ii=1:f.nblock, x = Gsb{ii}' * y; end, cpu etoc 'Gs time'
	cpu etic, for ii=1:f.nblock, x = G2b{ii}' * y; end, cpu etoc 'G2 time'
%profile on
%profile report
end
