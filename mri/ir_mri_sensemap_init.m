 function [sinit, bodycoil] = ir_mri_sensemap_init(ykj, varargin)
%function [sinit, bodycoil] = ir_mri_sensemap_init(ykj, varargin)
%|
%| Form initial sensitivity map estimates from surface coils and bodycoil.
%|
%| in
%| ykj [(N) ncoil] surface coil data (2D or 3D)
%|
%| option
%| 'type' initialization method (def: '' = 'ratio')
%|		options: 'median', 'order1', 'order2', ...
%|		if numeric, then simply return as "sinit"
%| 'bodycoil' [(N)] body coil data (def: SSoS of ykj with phase of 1st coil)
%| 'thresh' fraction of max |bodycoil| that is "good" (def: 0.05)
%| 'mask' [(N)] optional mask to pick good values instead of thresholding
%| 'chat' 0|1 verbosity
%|
%| out
%| sinit [(N) ncoil] initial sensitivity map estimates
%| bodycoil [(N)] computed SSoS bodycoil image (if not provided)
%|
%| used in ir_mri_sensemap_admm - see that m-file for details.
%|
%| Michael Allison
%| 2015-08-10 Jeff Fessler revise args, add tests, support 3D

if nargin < 1, ir_usage(), end
if streq(ykj, 'test', 4), ir_mri_sensemap_init_test(ykj), return, end

arg.type = '';
arg.thresh = 0.05;
arg.mask = [];
arg.bodycoil = [];
arg.chat = false;
arg = vararg_pair(arg, varargin);

bodycoil = arg.bodycoil;
if isempty(bodycoil)
	if ndims(ykj) < 3, fail('ykj must be 3D or 4D to infer bodycoil'), end
	ndim = ndims(ykj) - 1; % 2 or 3
	bodycoil = sqrt(sum(abs(ykj).^2, ndim+1)); % sum-of-squares
        tmp = stackpick(ykj, 1); % ykj(:,...,:,1)
        tmp = angle(tmp); % crude estimate of phase of image f_j
        bodycoil = bodycoil .* exp(1i * tmp);
else
	ndim = ndims(bodycoil);
	dim = size(bodycoil);
end

if ~isempty(arg.type) && isnumeric(arg.type) % trivial case
	sinit = arg.type;
return
end

dim = size(bodycoil);
ncoil = numel(ykj) / prod(dim);
if round(ncoil) ~= ncoil, fail('size'), end

if ~isempty(arg.mask)
	jf_equal(size(arg.mask), dim)
	if arg.chat
		warn('using mask instead of threshold to determine good pixels')
	end
	good = mask;
end

if arg.chat
	tell = @(s) printm(s);
else
	tell = @(s) nop(s);
end

sinit = zeros([prod(dim) ncoil], 'single');
for ic = 1:ncoil
	% cannot use stackpick because of case of just one coil
	if ndim == 2
		zj = ykj(:,:,ic);
	elseif ndim == 3
		zj = ykj(:,:,:,ic);
	end

	tmp = div0(zj, bodycoil); % usual ratio

	% determine "good" pixels
	if isempty(arg.mask)
		good = abs(bodycoil) > arg.thresh * max(abs(bodycoil(:)));
	end

	switch arg.type
	case {'', 'ratio'}
		% set all uncertain map values to median of good ones
		tell('Using zero background.')
		tmp(~good) = 0;

	case 'median'
		% set all uncertain map values to median of good ones
		tell('Using median background.')
		tmp1 = median(abs(tmp(good)));
		tmp2 = median(angle(tmp(good))); % dubious
		tmp(~good) = tmp1 .* exp(1i*tmp2);

	case 'avg'
		% set all uncertain map values to mean of good ones
		tell('Using mean background.')
		tmp1 = mean(abs(tmp(good)));
		tmp2 = mean(angle(tmp(good))); % dubious
		tmp(~good) = tmp1 .* exp(1i*tmp2);

	case 'zeros'
		tmp = zeros(dim);

	case 'order1'
		tell('Using 1st-order fit for background.')
		expo = [0 0; 1 0; 0 1];
		tmp = ortho_init(dim, expo, zj, bodycoil, good);

	case 'order2'
		tell('Using 2nd order fit for background.')
		expo = [0 0; 1 0; 0 1; 1 1; 2 0; 0 2];
		tmp = ortho_init(dim, expo, zj, bodycoil, good);

	case 'order3'
		tell('Using 3rd order fit for background.')
		expo = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1; 0 2; 1 2; ...
			3 0; 0 3];
		tmp = ortho_init(dim, expo, zj, bodycoil, good);

	case 'order4'
		tell('Using 4th order fit for background.')
		expo = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1; 0 2; 1 2; ...
			3 0; 0 3; 2 2; 3 1; 1 3; 4 0; 0 4];
		tmp = ortho_init(dim, expo, zj, bodycoil, good);

	otherwise
		fail('Unknown initialization type "%s"', init)
	end

	sinit(:,ic) = single(tmp(:));
end % ic
sinit = reshape(sinit, [dim ncoil]);

end %


function nop(varargin)
end % nop


% ortho_init()
% initialize using orthogonal polynomial basis fit
function init = ortho_init(dim, expo, zj, bodycoil, good)
switch numel(dim)
case 2
	%xx = ndgrid_jf('mat', 0:nx-1, 0:ny-1);
	A = cheby2d(dim(1), dim(2), expo);
case 3
	nx = dim(1);
	ny = dim(2);
	ny = dim(2);
	A = cheby3d(dim(1), dim(2), dim(3), expo);
otherwise
	fail('#dim = %d not done', numel(dim))
end

A = reshape(A, [prod(dim), size(expo,1)]);
maskBodyI = good .* bodycoil;
Aw = repmat(maskBodyI(:), [1 size(expo,1)]) .* A;
theta = pinv(Aw) * zj(:);
clear Aw
init = A * theta; % not Aw
init = reshape(init, size(bodycoil));

end % ortho_init()


% cheby2d()
function polys = cheby2d(nx,ny,expo)
nbases = length(expo);
x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
polys = zeros(nx,ny,nbases);
for n = 1:nbases
	tmp = polyval(ir_chebyshev_poly(expo(n,1)),x);
	tmp = tmp';
	X = repmat(tmp,[1 ny]);
	tmp = polyval(ir_chebyshev_poly(expo(n,2)),y);
	Y = repmat(tmp,[nx 1]);
	polys(:,:,n) = X .* Y;
end

end % cheby2d()


% cheby3d()
function polys = cheby3d(nx,ny,expo)
nbases = size(expo,1);
x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
z = linspace(-1,1,nz);
polys = zeros(nx, ny, nz, nbases, 'single');
for n = 1:nbases
	tmp = polyval(ir_chebyshev_poly(expo(n,1)), x);
	X = repmat(tmp', [1 ny]);
	tmp = polyval(ir_chebyshev_poly(expo(n,2)), y);
	Y = repmat(tmp, [nx 1]);
	tmp = polyval(ir_chebyshev_poly(expo(n,3)), z);
	Z = repmat(tmp, [nx 1]);
	polys(:,:,:,n) = X .* Y .* Z;
end

end % cheby3d()


% ir_mri_sensemap_init_test2()
% 2D
function ir_mri_sensemap_init_test2
nx = 28; ny = 32;
ig = image_geom('nx', nx, 'ny', ny, 'fov', 22);
tmp = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
%tmp = ones(nx,ny);
xtrue = tmp .* exp(1i * (tmp-2));
im(abs(xtrue))
ncoil = 4;
strue = ir_mri_sensemap_sim('nx', ig.nx, 'ny', ig.ny, 'dx', ig.dx, ...
	'ncoil', ncoil);
scale = 1 / sqrt(sum(abs(strue(end/2,end/2,:).^2)));
strue = strue * scale; % normalize so SSoS=1 near center
rng(0)
ykj = strue .* repmat(xtrue, [1 1 ncoil]) ...
	+ 2^-8 * (randn(size(strue)) + 1i * randn(size(strue)));
im plc 5 2
clim = minmax(abs(strue))';
im(1, 'row', 1, abs(strue), clim, 'true'), cbar
im(2, 'row', 1, abs(ykj)), cbar
sinit1 = ir_mri_sensemap_init(ykj);
im(3, 'row', 1, abs(sinit1), clim), cbar
slist = {'ratio', 'median', 'avg', 'order1', 'order2', 'order3', 'order4', 'zeros'};
for is=1:numel(slist)
	stype = slist{is};
	sinit2 = ir_mri_sensemap_init(ykj, 'type', stype);
	im(2+is, 'row', 1, abs(sinit2), clim, stype), cbar
	drawnow
end

end % ir_mri_sensemap_init_test2()


% ir_mri_sensemap_init_test3()
% 3D
function ir_mri_sensemap_init_test3
nx = 28; ny = 32; nz=6;
ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'fov', 22, 'zfov', 10);
tmp = ellipsoid_im(ig, 'shepp-logan-e3d', 'oversample', 2);
%tmp = ig.ones;
xtrue = tmp .* exp(1i * tmp); % arbitrary phase
%im(abs(xtrue))
nring = 3;
ncoil = 4 * nring;
strue = ir_mri_sensemap_sim('nx', ig.nx, 'ny', ig.ny, 'nz', ig.nz, ...
	'dx', ig.dx, 'nring', nring, 'ncoil', ncoil);
%im('row', 3, abs(strue))
scale = 1 / sqrt(sum(abs(strue(end/2,end/2,end/2,:).^2)));
strue = strue * scale; % normalize so SSoS=1 near center
rng(0)
ykj = strue .* repmat(xtrue, [1 1 1 ncoil]) ...
	+ 1 * 2^-12 * (randn(size(strue)) + 1i * randn(size(strue)));
im plc 4 1
clim = minmax(abs(strue))';
im(1, 'row', nring, abs(strue), clim, 'true'), cbar
titlef('True smaps for 3D [%d,%d,%d] with %d coils (%d rings of %d coils)', ...
	nx, ny, nz, ncoil, nring, ncoil/nring)
im(2, 'row', nring, abs(ykj), 'surface coils'), cbar
sinit1 = ir_mri_sensemap_init(ykj);
im(3, 'row', nring, abs(sinit1), clim, 'ratio'), cbar
slist = {'median'}; % 'avg', 'order1', 'order2', 'order3', 'order4', 'zeros'};
for is=1:numel(slist)
	stype = slist{is};
	sinit2 = ir_mri_sensemap_init(ykj, 'type', stype);
	im(3+is, 'row', nring, abs(sinit2), clim, stype), cbar
	drawnow
end
end % ir_mri_sensemap_init_test3()


% ir_mri_sensemap_init_test()
function ir_mri_sensemap_init_test(arg)
switch arg
case 'test2'
	ir_mri_sensemap_init_test2
case 'test3'
	ir_mri_sensemap_init_test3
case 'test'
	ir_mri_sensemap_init_test2
	if im, prompt, end
	ir_mri_sensemap_init_test3
otherwise
	fail('unknown test "%s"', arg)
end
end % ir_mri_sensemap_init_test()
