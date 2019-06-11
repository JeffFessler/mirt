 function ugibb = ugibb_form(pic, varargin)
%function ugibb = ugibb_form(pic, [options])
%|
%| form a "ugibb" structure from a 2D or 3D picture,
%| of potential use for regularization with side information.
%| Specifically can be used with ASPIRE 'ugibb' regularization option.
%| See Fessler et al, IEEE Tr. Nuc. Science 39(5):1464, Oct. 1992,
%| "Regularized emission image reconstruction using imperfect side information"
%| and Comtat et al, PMB 47(1):1-20, Jan. 2002, "Clinically feasible
%| reconstruction of 3d whole-body PET/CT data using blurred anatomical labels"
%|
%| in
%|	pic	[(N)]		2D or 3D array of region labels
%|				usually 0=background, 1=roi1, 2=roi2, ...
%|
%| options
%|	'threshold'		default: 0 (preserve any differences)
%|	'offsets'		cf penalty_offsets.m
%|	'blur_fwhm'		for blurred masks (ala Comtat); pic must be 4D
%|
%| out
%|	ugibb	[(N) 4 or 5]	3D [nx ny 4] or 4D [nx ny nz 5] weights
%|
%| Copyright 2004-11-14, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(pic, 'test'), ugibb_test, return, end

arg.threshold = 0;
arg.offsets = [];
arg.blur_fwhm = 0; % for blurred labels [pixels]
arg = vararg_pair(arg, varargin);

switch ndims(pic)
case 2
	ugibb = ugibb_form_2d(pic, arg);
case 3
	ugibb = ugibb_form_3d(pic, arg);
case 4
	ugibb = ugibb_form_3d_blur(pic, arg);
otherwise
	fail('unknown')
end


% ugibb_form_2d()
% 2d case
%
function ugibb = ugibb_form_2d(pic, arg)
[nx ny] = size(pic);
if isempty(arg.offsets), arg.offsets = [1 nx nx+1 nx-1]; end
ii = (nx+2):(nx*ny); % trick: no dx on 1st image row...

if 1
	dd = zeros(nx*ny, length(arg.offsets));
	for io=1:length(arg.offsets)
		off = arg.offsets(io);
%		dd(:,ii) = pic(ii) - pic(ii-off);
		ii = (1+off):(nx*ny);
		dd(ii,io) = pic(ii) - pic(ii-off);
	end
	dd = reshape(dd, nx, ny, length(arg.offsets));
else
	dx = zeros(nx,ny); dy = dx; dp = dx; dn = dx;
	dx(ii) = pic(ii) - pic(ii-1);
	dy(ii) = pic(ii) - pic(ii-nx);
	dp(ii) = pic(ii) - pic(ii-nx+1);
	dn(ii) = pic(ii) - pic(ii-nx-1);
	dd(:,:,1) = dx;
	dd(:,:,2) = dy;
	dd(:,:,3) = dp;
	dd(:,:,4) = dn;

end
ugibb = single(abs(dd) < arg.threshold);

% fix the borders (where no neighbors exist)
if 1
	ugibb(1,:,[1 4]) = 0;
	ugibb(:,1,[2 3 4]) = 0;
	ugibb(end,:,3) = 0;
end

% ugibb(:,:,[3 4]) = ugibb(:,:,[3 4]) / sqrt(2);


% ugibb_form_3d()
% 3d case
% with careful treatment of edges so zero weight there!
%
function ugibb = ugibb_form_3d(pic, arg)
[nx ny nz] = size(pic);

if isempty(arg.offsets)
	arg.offsets = [1 nx nx+1 nx-1 nx*ny];
end
dx = inf(nx,ny,nz,'single'); dy = dx; dp = dx; dn = dx; dz = dx;

ix = col([2:nx]'*ones(1,ny*nz) + ones(nx-1,1)*([1:ny*nz]-1)*nx);
dx(ix) = pic(ix) - pic(ix-1);

[ix iy iz] = ndgrid(1:nx,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
dy(ii) = pic(ii) - pic(ii-nx);

[ix iy iz] = ndgrid(1:nx-1,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
dp(ii) = pic(ii) - pic(ii-nx+1);

[ix iy iz] = ndgrid(2:nx,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
dn(ii) = pic(ii) - pic(ii-nx-1);

iz = (nx*ny+1):(nx*ny*nz);
dz(iz) = pic(iz) - pic(iz-nx*ny);

if 0 % old way
	dd = zeros(nx,ny,5,nz);
	dd(:,:,1,:) = dx;
	dd(:,:,2,:) = dy;
	dd(:,:,3,:) = dp;
	dd(:,:,4,:) = dn;
	dd(:,:,5,:) = dz;
else
	dd = cat(4, dx, dy, dp, dn, dz); % [(N) 5]
end

ugibb = abs(dd) < arg.threshold; % logical
ugibb = single(ugibb); % caution: needed for diagonal scaling to work! 2009-9-22

% scale diagonals
if 0 % old
	ugibb(:,:,[3 4],:) = ugibb(:,:,[3 4],:) / sqrt(2);
else
	ugibb(:,:,:,[3 4]) = ugibb(:,:,:,[3 4]) / sqrt(2);
end

if 0
	rot180 = @(x) rot90(x, 2);
	dx = flipud(conv2(flipud(pic), [1 -1]', 'same'));
	dy = fliplr(conv2(fliplr(pic), [1 -1], 'same'));
	dn = rot180(conv2(rot180(pic), [1 0; 0 -1], 'same'));
	dp = fliplr(conv2(fliplr(pic), [1 0; 0 -1], 'same'));

	da = dx;	% all 4 of them
	da(:,:,2) = dy;
	da(:,:,3) = dp;
	da(:,:,4) = dn;
	u = abs(da) < arg.threshold;
end


% ugibb_form_3d_blur()
% 3d case ala Comtat with optional blur
% with careful treatment of edges so zero weight there!
% masks should be [nx ny nz nr] where "nr" is number of regions
%
function ugibb = ugibb_form_3d_blur(masks, arg)
[nx ny nz nr] = size(masks);

if isempty(arg.offsets)
	arg.offsets = [1 nx nx+1 nx-1 nx*ny];
end

if arg.threshold ~= 0
	fail 'threshold not relevant to blurred label approach'
end

if arg.blur_fwhm ~= 0
	kern = gaussian_kernel(arg.blur_fwhm);
	if arg.blur_fwhm ~= 0
		masks = single(masks); % no longer logical
	end
	nk = length(kern);
	for ir=1:nr
		tmp = masks(:,:,:,ir);
		tmp0 = tmp(1); % border value
		tmp = tmp - tmp0; % try to make border 0
		tmp = convn(tmp, reshape(kern, [nk 1 1]), 'same');
		tmp = convn(tmp, reshape(kern, [1 nk 1]), 'same');
		tmp = convn(tmp, reshape(kern, [1 1 nk]), 'same');
		tmp = tmp + tmp0; % restore original DC level
		masks(:,:,:,ir) = tmp;
	end
end

dx = zeros(nx,ny,nz,'single'); dy = dx; dp = dx; dn = dx; dz = dx;

ix = col([2:nx]'*ones(1,ny*nz) + ones(nx-1,1)*([1:ny*nz]-1)*nx);
sum = 0;
for ir=1:nr
	tmp = masks(:,:,:,ir);
	sum = sum + tmp(ix) .* tmp(ix-1);
%	dx(ix) = pic(ix) - pic(ix-1);
end
dx(ix) = sum;

[ix iy iz] = ndgrid(1:nx,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
sum = 0;
for ir=1:nr
	tmp = masks(:,:,:,ir);
	sum = sum + tmp(ii) .* tmp(ii-nx);
%	dy(ii) = pic(ii) - pic(ii-nx);
end
dy(ii) = sum;

[ix iy iz] = ndgrid(1:nx-1,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
sum = 0;
for ir=1:nr
	tmp = masks(:,:,:,ir);
	sum = sum + tmp(ii) .* tmp(ii-nx+1);
%	dp(ii) = pic(ii) - pic(ii-nx+1);
end
dp(ii) = sum;

[ix iy iz] = ndgrid(2:nx,2:ny,1:nz);
ii = ix + ((iz-1)*ny+(iy-1))*nx;
sum = 0;
for ir=1:nr
	tmp = masks(:,:,:,ir);
	sum = sum + tmp(ii) .* tmp(ii-nx-1);
%	dn(ii) = pic(ii) - pic(ii-nx-1);
end
dn(ii) = sum;

iz = (nx*ny+1):(nx*ny*nz);
sum = 0;
for ir=1:nr
	tmp = masks(:,:,:,ir);
	sum = sum + tmp(iz) .* tmp(iz-nx*ny);
%	dz(iz) = pic(iz) - pic(iz-nx*ny);
end
dz(iz) = sum;

ugibb = cat(4, dx, dy, dp, dn, dz); % [(N) 5]

% scale diagonals
ugibb(:,:,:,[3 4]) = ugibb(:,:,:,[3 4]) / sqrt(2);


%
% ugibb_test.m
%
% test ugibb_form() and ugibb_sum()
function ugibb_test

if 0 % 2d
	pic = zeros(7,5);
	pic(3,2) = 1;
	pic = zeros(8,10);
	pic(4:6,4:6) = 1;
	[nx,ny] = size(pic);
else % 3d
	pic = zeros(12,10,6);
	pic(3:6,3:6,2:4) = 1;
	pic(4:5,4:5,3) = 2;
	pic(end/2,end/2,end) = 1;
end

if 0
	u = ugibb_form(pic, 'threshold', 0.02);
else
	masks = cat(4, pic==0, pic==1, pic==2);
%	clf, im row 3, im(masks), return
	u = ugibb_form(masks, 'blur_fwhm', 0.1);
end
s = ugibb_sum(u);

im plc 2 4
p = @(ii) stackpick(u, ii);
im(2, p(1), 'wx'), cbar h
im(3, p(2), 'wy'), cbar h
im(6, p(3), 'wp'), cbar h
im(7, p(4), 'wn'), cbar h
im(8, p(5), 'wz'), cbar h
im(1, pic, 'pic'), cbar h
im(5, s, 'sum'), cbar h
