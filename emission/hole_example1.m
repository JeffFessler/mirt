% hole_example1.m
%
% 2D example showing the "hole in sphere" effect for OSEM SPECT
% and various regularization methods that may help overcome it.
%
% Copyright 2004-5-13, Jeff Fessler, University of Michigan

if ~has_aspire, return, end

if ~isvar('A0'), printm 'generate system matrix'
	f.dir = test_dir;
	f.dsc0 = [f.dir 't0.dsc'];
	f.wtf0 = strrep(f.dsc0, 'dsc', 'wtf');
	if 1 % todo: replace with aspire_pair
		os_run(['wt -chat 0 dsc 12 fwhm_detector 5 >! ' f.dsc0])
		os_run(['echo y | wt -chat 0 gen ' f.dsc0 ' row'])
	end

	A0 = Gtomo2_wtmex(f.wtf0);
	ig = A0.arg.ig;
	sg = A0.arg.sg;
	nx = ig.nx; nb = sg.nb;
	ny = ig.ny; na = sg.na;
end

if ~isvar('Ab0'), printm 'Ab'
	f.nblock = 6;
	Ab0 = Gblock(A0, f.nblock);
	Ab20 = Gblock(A0, 20);
end


if ~isvar('yi'), printm 'projection data'
%	xtrue = [6.5 -0.5 5 5 0 100]; over = {'oversample', 4};
%	xtrue = ellipse_im(ig, xtrue, over{:});
	over = {'oversample', 1};
	xtrue1 = ellipse_im(ig, [18.5 -0.5 5 5 0 100], over{:});
	xtrue2 = ellipse_im(ig, [-3.5 -0.5 3 3 0 100], over{:});
	xtrue3 = ellipse_im(ig, [-18.5 -0.5 2 2 0 100], over{:});
	xtrue = xtrue1 + xtrue2 + xtrue3;
	im(xtrue, 'xtrue'), cbar
prompt

	ri = 1;
	yi = A0 * xtrue + ri;
	im(yi, 'yi'), cbar
prompt
end


% uniform initial image
xinit = ig.ones;

if 0 && ~isvar('xfbp'), printm 'FBP'
	tmp = fbp2(sg, ig);
	xfbp = fbp2(yi-ri, tmp, 'window', 'hann');
	im(xfbp, 'fbp'), cbar
prompt
end


if ~isvar('xe0'), printm 'EM'
	f.niter_em = 5;
%	xe0 = eml_em(xinit(ig.mask), A0, yi(:), 1, ri, [], f.niter_em);
	xe0 = eml_em(xinit(ig.mask), A0, yi(:), 1, ri, ...
		'isave', 'all', 'niter', f.niter_em);
	xe0 = ig.embed(xe0);
	im(xe0, 'EM'), cbar
prompt
end

xinit = xe0(:,:,end); % last EM iter


% OS-EM iterations
if ~isvar('xo0'), printm 'run os-em'
	f.niter = 10; % for testing
	if isempty(caller_name)
		f.niter = 200; % for full comparison
	end

	rri = ri * ones(nb,na);
	xo0 = eml_osem(xinit(ig.mask), Ab0, yi, [], rri, f.niter);
	xo0 = ig.embed(xo0);
	im(xo0(:,:,end), 'OSEM'), cbar
prompt
end


% regularization (edge preserving)
if ~isvar('R'), printm 'R'
	f.l2b = -9;
	f.l2b = -11; % for binary case
	f.potential = 'hyper3'; f.delta = 3;
%	f.potential = 'cauchy'; f.delta = 10; % this one used for grant
	f.potential = 'cauchy'; f.delta = 2;
%	R = Robject(ig.mask, 'type_denom', 'matlab', ...
%		'beta', 2^f.l2b, 'potential', f.potential, 'delta', f.delta);
	R = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'beta', 2^f.l2b, 'pot_arg', {f.potential, f.delta});
	clear xp0
prompt
end

% Incremental EM-3, aka C-OS-3
if ~isvar('xp0'), printm 'C-OS-3'
%	xp0 = epl_os_emdp(xinit(ig.mask), Ab0, yi, ones(size(yi)), rri, R, ...
%		0+1*f.niter, 200, 1);
	xp0 = epl_inc(xinit(ig.mask), Ab20, yi, ones(size(yi)), rri, R, ...
		'niter', f.niter, 'hds', 3, 'os', 40, 'pixmax', 200, 'chat', 0);
	xp0 = ig.embed(xp0);
	im(xp0(:,:,end), 'C-OS-3')
prompt
end

% regularization with perfect side info
if ~isvar('Rs'), printm 'Rs'
	ug = ugibb_form(xtrue > 50, 'threshold', 0.5, 'offsets', R.offsets);
	im(ug)
	f.l2b = -9;
%	Rs = Robject(ig.mask, 'type_denom', 'matlab', ...
%		'beta', 2^f.l2b, 'user_wt', ug);
	Rs = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'beta', 2^f.l2b, 'user_wt', ug);
	clear xs0
prompt
end

if ~isvar('xs0'), printm 'recon with side info'
	xs0 = epl_inc(xinit(ig.mask), Ab20, yi, ones(size(yi)), rri, Rs, ...
		'niter', f.niter, 'hds', 3, 'os', 40, 'pixmax', 200, 'chat', 0);
	xs0 = ig.embed(xs0);
	im(xs0(:,:,end), 'xs0')
prompt
end


if 1 && im, printm 'show "hole in sphere"'
	iy = ny/2;
	t = squeeze(xo0(:,iy,:));
	t = squeeze(xp0(:,iy,:));
	im(t)

	ix = 1:nx;
%	ix = 10:60;
	clf
%	subplot(211)
%	plot(ix, xtrue(ix,iy), 'y-o', 1:nx, t);

	plot(	ix, xtrue(ix,iy), 'y-x', ...
		ix, xo0(ix,iy,end), 'c-o', ...
		ix, xp0(ix,iy,end), 'g.-', ...
		ix, xs0(ix,iy,end), 'r.-')
	legend('True', 'ML: OS-EM', 'PL: C-OS-3', 'PL: side', 4)
	xlabel 'Horizontal Pixel Index'
	ylabel 'Activity [arbitrary units]'
	title 'Profile'
	axisy([0 160])
prompt
end

if im, printm 'show total activity per iteration'
	t = xtrue3 > 0;
	sum_true = sum(xtrue(t));
	t1 = reshape(xo0, [], f.niter+1); t1 = sum(t1(t,:)) / sum_true;
	t2 = reshape(xp0, [], f.niter+1); t2 = sum(t2(t,:)) / sum_true;
	t3 = reshape(xs0, [], f.niter+1); t3 = sum(t3(t,:)) / sum_true;
	ii = 0:f.niter;
	clf, % subplot(212)
	plot( ...
		ii(1), t3(1), 'y-s', ...
		ii(1), t2(1), 'g.-', ...
		ii(1), t1(1), 'c--o', ...
		ii(1:8:end), t3(1:8:end), 'ys', ...
		ii(1:8:end), t2(1:8:end), 'g.', ...
		ii(1:8:end), t1(1:8:end), 'co', ...
		ii, t3, 'y-', ...
		ii, t2, 'g-', ...
		ii, t1, 'c--', ...
		ii, 1, 'y:')
	legend('PL side', 'PL: C-OS-3', 'ML: OS-EM', 4)
%	legend('PL: C-OS-3 with side information', 'ML: OS-EM', 4)
	xlabel 'Iteration', ylabel 'Activity Recovery'
	axis([0 120 0.8 1.02])
	axis([0 120 0.6 1.02])
%	ir_savefig hole_example_recovery1 c
%	set(gcf, 'InvertHardCopy', 'off'); print -djpeg90 hole_example_recovery2
prompt
end

if 1 && im, printm 'examine "hole in sphere"'
	t = xtrue1(:) > 99;
	t1 = reshape(xo0, [nx*ny f.niter+1]);
	t1 = t1(t, :);
	t2 = reshape(xp0, [nx*ny f.niter+1]);
	t2 = t2(t, :);
	ii = 0:f.niter;
	clf, % subplot(212)
	plot(	...
		ii(1), max(t1(:,1), [], 1), 'c-o', ...
		ii(1), min(t1(:,1), [], 1), 'c--x', ...
		ii(1), max(t2(:,1), [], 1), 'g-', ...
		ii(1), min(t2(:,1), [], 1), 'g--', ...
		ii(1:8:end), max(t1(:,1:8:end), [], 1), 'co', ...
		ii(1:8:end), min(t1(:,1:8:end), [], 1), 'cx', ...
...%		ii(1:8:end), max(t2(:,1:8:end), [], 1), 'g-', ...
...%		ii(1:8:end), min(t2(:,1:8:end), [], 1), 'g--', ...
		ii, max(t1, [], 1), 'c-', ...
		ii, min(t1, [], 1), 'c--', ...
		ii, max(t2, [], 1), 'g-', ...
		ii, min(t2, [], 1), 'g--', ...
		ii, 100, 'y:')
	legend('ML: OSEM max', 'ML: OSEM min', 'PL: C-OS-3 max', 'PL: C-OS-3 min')
	xlabel 'Iteration'
	ylabel 'Activity [arbitrary units]'
	title 'Maximum and Minimum Values within ROI'
	axisy([0 165])
prompt
end

if 1, printm 'show reconstructions' % for Ken
	ix = 1:nx; iy=1:ny;
%	ix = 7+[nx/4:(3*nx/4)]; iy=ny/4:(3*ny/4);
	t = stackup(xtrue(ix,iy), xo0(ix,iy,end), xp0(ix,iy,end), xs0(ix,iy,end));
	im clf, im(t, [0 120]), cbar
%	t(t>120) = 120; t = uint8(255 * t/120);
%	imwrite(t', 'hole.tif')
end

% ir_savefig fig_hole_example1_new
