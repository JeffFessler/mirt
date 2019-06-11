  function [C, rj] = penalty2_nuyts(type, wang, ang, mask)
%|function [C, rj] = penalty2_nuyts(type, wang, ang, mask)
%|
%| Design the penalty matrix "C" (and penalty coefficients "rj")
%| for a quadratic penalty R(x) = 1/2 x' C * C * x with:
%| 1st-order differences and a 2nd-order neighborhood (8 neighbors).
%| Design based on heuristics in Nuyts and Fessler, IEEE T-MI, 2003.
%| Actually, it is based on Nuyts' IDL code.
%| For 2D parallel-beam tomography.
%|
%| in
%|	type			'test' or 'quad,d1,n2'
%|	wang	[na nx ny]	angular weights for each pixel
%|	ang	[na]		projection angles (in radians)
%|	mask	[nx ny]		(logical) reconstruction support
%|
%| out
%|	rj	[nx ny 4]	horiz,vert,diag1,diag2 weights
%|
%| caution: sqrt(2) factors are already built into C (via rj).
%|
%| Copyright 2003-5-23, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

% self-test compares Nuyts design to Fourier-NNLS design.
if nargin == 1 && streq(type, 'test')
	penalty2_nuyts_test
return
end

%arg.edge_type = 'leak'; % this could be a user-selectable option.

switch type
case 'quad,d1,n2'
	if nargin ~= 4, error 'usage', end
	[C rj] = penalty2_nuyts_quad_d1_n2(type, wang, ang, mask);
otherwise
	error(['unknown type: ' type])
end


% penalty2_nuyts_quad_d1_n2()
% Nuyts design of 1st-order difference, 2nd-order neighborhood
function [C, rj] = penalty2_nuyts_quad_d1_n2(type, wang, ang, mask)

[na nx ny] = size(wang);
np = nx*ny;
wang = reshape(wang, [na np]);
if length(ang) ~= na, error 'ang size', end

if ~isvar('mask') || isempty(mask)
	mask = true(nx,ny);
else
	if any(size(mask) ~= [nx ny]), error 'mask dims', end
end

% match each angle to one of the 4 angular blocks
ia = round(ang / (pi/4));
ia = 1 + mod(ia, 4); % {1,2,3,4}

% sum the "information" values (certainties) in each block
order = [1 3 4 2]; % 0 90 45 135 ... fix: fails if flipy=-1?
fish = zeros(4,np); % oriented fisher info for each pixel
for ii=1:4
	fish(ii,:) = sum(wang(ia==order(ii),:),1) / na;
end

if 0 % compare to Johan's results!
	joh = fld_read('~fessler/ocpl/fisher_theta_before.fld');
	joh = joh(:,:,[1 3 4 2]);
	joh(:,:,1) = fliplr(joh(:,:,1));
	joh(:,:,2) = fliplr(joh(:,:,2));
	joh(:,:,3) = fliplr(joh(:,:,3));
	joh(:,:,4) = fliplr(joh(:,:,4));

	tmp = reshape(fish, [4 nx ny]);
	tmp = permute(tmp, [2 3 1]);
	p = 0.2;
	clf, im(121, 'colorneg', tmp .^p), im(122, 'colorneg', joh .^p)

	if 0 % use johan's fisher_theta ?
		tmp = permute(joh, [3 1 2]);
		fish = reshape(tmp, [4 np]);
	end
end

% nuyts heuristic modifications
fmax = max(fish, [], 1); % [1 np]
hv = (fmax == max(fish([1 2],:), [], 1)); % horiz./vert. dominate

% where horizontal/vertical dominate, no diagonal components
rj = zeros(4,np);
rj([1 2],hv) = fish([1 2],hv);

% where diagonals dominate, include some horizontal/vertical smoothing
fmin = min(fish([3 4],~hv), [], 1); % min of diagonals
rj([1 2], ~hv) = [fmin; fmin];
rj([3 4], ~hv) = fish([3 4],~hv) - [fmin; fmin];

% 2d gaussian blur
kern = gaussian_kernel(3); % FWHM=3 per Nuyts code
rj = permute(reshape(rj, [4 nx ny]), [2 3 1]); % [nx ny 4]
for ii=1:4
	rj(:,:,ii) = conv2(rj(:,:,ii), kern, 'same');
	rj(:,:,ii) = conv2(rj(:,:,ii), kern', 'same');
end

% normalize to preserve total fisher information
fish_total = reshape(sum(fish,1), [nx ny]);
rj = 2 * rj .* repeat_slice(fish_total ./ (sum(rj, 3) + eps), 4);

% nuyts divides diagonal coefficients by sqrt(2), as many others
% have done, although Fourier analysis would question this practice...
rj(:,:,[3 4]) = rj(:,:,[3 4]) / sqrt(2);

C = C2sparse('leak', ones(nx,ny), 8, 0, 0);
i1 = 1:(nx*ny);
C = [ ...
	spdiag(sqrt(col(rj(:,:,1))), 'nowarn') * C(0*nx*ny+i1,:);
	spdiag(sqrt(col(rj(:,:,2))), 'nowarn') * C(1*nx*ny+i1,:);
	spdiag(sqrt(col(rj(:,:,3))), 'nowarn') * C(2*nx*ny+i1,:);
	spdiag(sqrt(col(rj(:,:,4))), 'nowarn') * C(3*nx*ny+i1,:)];


% penalty2_nuyts_test()
% test and compare. this is longer than the routine to be tested...
function penalty2_nuyts_test

ig = image_geom('nx', 100, 'ny', 98, 'dx', 4.0);
% from Nuyts paper:
sg = sino_geom('par', 'nb', 102, 'na', 80, 'dr', ig.dx, 'strip_width', 'd');
A = Gtomo2_wtmex(sg, ig, 'grouped', 'col'); % system matrix

a = 8 / ig.dx;
xtrue = ellipse_im(ig, [0 0	45 45 0 a/4] * ig.dx, 'oversample', 4) + ...
	ellipse_im(ig, ...
		[-15 15	14 14	0 a;
		-15 15	12 12	0 -a;
		25 0	14 14	0 a;
		25 0	12 12	0 -a] * ig.dx); % emission object

tmp = ellipse_im(ig, [[ig.dx ig.dy]*0.5 [1 1]*4.5*ig.dx 0 1], 'oversample', 1);
tmp = tmp(end/2+1+[-4:4],end/2+1+[-4:4]) > 0; % [9 9]
tmp = tmp / sum(tmp(:));
misk = conv2(single(xtrue > 0), single(tmp), 'same') >= 1-1e-5; % interior

% (nonuniform) sinogram weights w_i
xmu = 0.0095 * (xtrue ~= 0); % uniform attenuation
ci = exp(-(A * xmu));
yt = ci .* (A * xtrue)/4 + 1; % include additive background
wi = ci.^2 ./ max(yt,1);
W = diag_sp(wi(:));
clear xmu ci yt

% stayman angular weighting factors for each pixel
wang = A.arg.stayman2_factors(A, wi); % [na nx ny]

kappa = sqrt(squeeze(mean(wang,1))); % [nx ny] certainties

cx = ig.nx/2+1;
cy = ig.ny/2+1;
beta = 2^14; % 3.8 FWHM

%offset = [35 0];
%offset = [-20 -20]; % wiggly 4
%offset = [-15 -15]; % nice
%offset = [0 0];
%offset = [25 0]; % right circle center
%offset = [23 23];
offset = [20 10]; % dramatic!

do_psf = @(C, beta, W) ...
	qpwls_psf(A, C, beta, ig.mask, W, 'offset', offset, 'shift0', 1);

% target PSF for QPULS (unweighted)
C0 = C2sparse('leak', ig.mask, 4, 0);
psf.target = do_psf(C0, beta, 1);
%clf, im(psf.target), cbar, return

% center scaling only regularization (standard)
C = double(kappa(cx,cy)) * C0;
psf.center = do_psf(C, beta, W);

% certainty-weighted (Fessler & Rogers, 1996)
[C tmp] = C2sparse('leak', kappa, 4, 0);
C = spdiag(sqrt(tmp), 'nowarn') * C;
psf.cert = do_psf(C, beta, W);

% Nuyts & Fessler, 2003
% [nx,ny,4]
[C rj_nuyts] = penalty2_nuyts('quad,d1,n2', wang, sg.ar, ig.mask);
psf.nuyts = do_psf(C, beta, W);

% Fessler Fourier-NNLS method
[C rj_fess] = penalty2_design('quad,d1,n2', wang, sg.ar, ig.mask);
psf.fess = do_psf(C, beta, W);

if im
	clf, im pl 4 1
	tmp = permute(rj_fess .* repeat_slice(misk,4), [1 3 2]);
	im(1, reshape(tmp, [ig.nx*4, ig.ny]), 'Fourier-NNLS'), cbar
	hold on, plot(cx+offset(1)+[0:3]*ig.nx, cy+offset(2), 'y+'), hold off
	tmp = permute(rj_nuyts .* repeat_slice(misk,4), [1 3 2]);
	im(2, reshape(tmp, [ig.nx*4, ig.ny]), 'Nuyts'), cbar
	hold on, plot(cx+offset(1)+[0:3]*ig.nx, cy+offset(2), 'c+'), hold off

	jx = [-8:8] + cx; % don't need offset!
	jy = [-8:8] + cy; % for psf display

	lirs = [psf.target(jx,jy); psf.center(jx,jy);
		psf.cert(jx,jy); psf.nuyts(jx,jy);
		psf.fess(jx,jy)];
	type = {'Target', 'Standard', 'Certainty', 'Nuyts', 'Fourier-NNLS'};
	name = {'target', 'center', 'cert', 'nuyts', 'fess'};

	jjx = 1:length(type)*length(jx);
	jjy = cy-jy;
	ax = [minmax(jjx); minmax(jjy)];
	im(3, jjx, jjy, lirs), axis(ax), axis xy, axis equal, xtick off
	cbar
	title(sprintf('Local impulse response functions at (%d,%d)', ...
		offset(1), offset(2)))
	im subplot 4
	c = [0.1 0.25 0.5 0.75 0.9] * max(psf.target(:));
	contour(jjx, jjy, lirs', c, '-'), axis(ax), axis equal, xtick off
	title('Contours and FHWM')
	for ii=1:length(type)
		ix = 9+(ii-1)*17;
		text(ix, 6, type{ii}, 'horiz', 'center')
		text(ix, -6, sprintf('%4.2f', fwhm2(psf.(name{ii}))), 'horiz', 'center')
	end
end
