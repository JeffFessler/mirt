% Glsi_blur_test.m
% Test the Glsi_blur object

if ~isvar('G'),	disp 'setup Glsi_blur_test'
	psf = [0 1 2 1 0; 1 2 3 2 1; 0 1 2 1 0];
	psf = psf / sum(psf(:));
	nx = 128;
	ny = 100;

	G = Glsi_blur(nx, ny, psf, [], 1);
	nb = G.nb;
	na = G.na;
	clf, im(331, G.psf, 'psf')
	im(332, G.mask, 'mask')
prompt
end

% test G and G'
if 1
	x = shepplogan(nx, ny, 1);
	y1 = reshape(G * x(:), nb, na);
	y2 = reshape(G' * x(:), nb, na);
	im(334, y1, 'forward G')
	im(335, y2, 'forward A')
	im(336, y2-y1, 'A-G')
prompt
end

% make a noisy blurry image
if ~isvar('y')
	rng(0)
	x = shepplogan(nx, ny, 1);
	im(334, x, 'true image x')
	y = G * x;
	y = y + 0.5 * randn(size(y));	% add noise
	im(335, max(y,0), 'noisy blurred image y')
prompt
end

%
% penalty function
%
if ~isvar('R')
	f.l2b = 14;
	f.delta = 0.1;
	R = Robject(G.mask, 'potential', 'huber', ...
		'beta', 2^f.l2b, 'delta', f.delta, 'type_denom', 'matlab');
prompt
end

%
% apply iterative restoration
%

if 1
	xs = pwls_sps_os(max(y(:),0), y(:), 1, G, R, 15, 1e9);
	xs = embed(xs, G.mask);
	clf, im(xs)
prompt
end


if 0
	y = ones(nb, na);
	x1 = reshape(G' * y(:), nx, ny);
	x2 = reshape(A' * y(:), nx, ny);
	im(337, x1, 'backward G')
	im(338, x2, 'backward A')
	im(339, x2-x1, 'A-G')
end
