function [xs, xc] = sense_recon(ysense, smap, varargin)

% input R?
% warn user if their R means throwing away data?
% ysense = fully samples with zeros?  
% allow 2D y w/ kspace input, as well as 4D y
% Note: this is for undersampling in 1 dim only (for now)


if (nargin == 1) && isequal(ysense, 'test')
    sense_recon_test();
    return;
end

if (nargin == 1) && isequal(ysense, 'test2')
    sense_recon_test2();
    return;
end

[ny, nz, nc] = size(smap);
[nky, nkz, nc2, M] = size(ysense);

dky = ny/nky;
dkz = nz/nkz;    % not used in current version, nz and nkz should always ==

% *** tmp ****
% dky = 1;

if nc ~= nc2
    error('data array must have same number of coils as smap');
end

% establish defaults
arg.R = 2;
arg.kspace = [];
arg.phi = 1;
arg = vararg_pair(arg, varargin);

% keyboard

xs = zeros(ny, nz, M);
for ifr = 1:M
    yfr = ysense(:,:,:, ifr);
    
    xc = zeros(nky, nz, nc);
    for ic = 1:nc
        % create aliased coil images 
        xc(:,:, ic) = arg.R * 1/dky * fftshift(ifft2(fftshift(yfr(:,:,ic))));
%         xc(:,:, ic) = fftshift(ifft2(fftshift(yfr(:,:,ic))));

    end
    
% disp('sense_recon: before call to desense, check xc size')
% keyboard

    xfr = desense(xc, smap, 'phi', arg.phi);

    xs(:,:, ifr) = xfr;
end

% disp('sense_recon: desense complete,check out coil images (xc) and such')
% keyboard

end


function sense_recon_test()

s = load('~/research/dce/PEsampling/savedData/psf_hs12_3.mat');
load('~/research/dce/PEsampling/savedData/whatTruirNeeds.mat', 'smap');

[ny, nz, nc] = size(smap);
M = 3;

ys_sm = s.yhs_psf(:, :, :, 1:M);
ys = zeros(ny, nz, nc, M);
ys(1:2:end, :, :, :) = ys_sm;

xs = sense_recon(ys, smap);
xs_sm = sense_recon(ys_sm, smap);


figure(1)
im('row', M, 211, xs_sm), cbar
im('row', M, 212, xs), cbar

% keyboard


end


function sense_recon_test2()

% s = load('~/research/dce/PEsampling/savedData/psf_hs12_3.mat');
%load('~/research/dce/PEsampling/savedData/whatTruirNeeds.mat', 'smap');
load('whatTruirNeeds.mat', 'smap');

[ny, nz, nc] = size(smap);
M = 2;

% 2 coils only
nc2 = 2;
smap2 = zeros(ny, nz, nc2);
ph2 = (-nz/2:nz/2-1)/nz * 2*pi;
ph2 = repmat(ph2, [ny 1]);
ph2 = zeros(ny, nz); % jf
%smap2(1:end/2, : , 1) = 1 .* exp(1i * ph2(1:end/2, :));
%%smap2(1:end/2, : , 1) = repmat(1 + [1:ny/2]'/(ny/2), [1 nz]);
smap2(1:end/2, : , 1) = (1 + [1:ny/2]'/(ny/2)) * (2 + [1:nz]/nz);
smap2(end/2+1:end, :, 2) = 1 .* exp(1i * ph2(end/2+1:end, :));

% 4 coils
nc4 = 4;
smap4 = zeros(ny, nz, nc4);
smap4( 1:end/4, :, 1) = 1;

smap4( 1*end/4+1:2*end/4, :, 2) = 0.5;
smap4( 3*end/4+1:end, :, 2) = 0.5;

smap4( 2*end/4+1:3*end/4, :, 3) = 0.5;
smap4( 1*end/4+1:2*end/4, :, 3) = 0.5;

smap4( 3*end/4+1:end, :, 4) = 0.5;
smap4( 2*end/4+1:3*end/4, :, 4) = 0.5;

figure(1), clf
im plc 3 2
im(1, 'row', M, smap), cbar
im(3, 'row', M, smap2), cbar
im(5, 'row', M, smap4), cbar
im(2, 'row', M, angle(smap)), cbar
im(4, 'row', M, angle(smap2)), cbar
im(6, 'row', M, angle(smap4)), cbar

ys = impulsedata(smap, M);
ys2 = impulsedata(smap2, M);
ys4 = impulsedata(smap4, M);

ys_sm = ys(1:2:end, :, :, :);
ys2_sm = ys2(1:2:end, :, :, :);
ys4_sm = ys4(1:2:end, :, :, :);

figure(2), clf
im plc 3 2
im(1, 'row', M, ys), cbar
im(3, 'row', M, ys2), cbar
im(5, 'row', M, ys4), cbar
im(2, 'row', M, ys_sm), cbar
im(4, 'row', M, ys2_sm), cbar
im(6, 'row', M, ys4_sm), cbar


[xs_sm xc_sm] = sense_recon(ys_sm, smap);
keyboard; return

[xs2_sm xc2_sm] = sense_recon(ys2_sm, smap2);
xs2 = sense_recon(ys2, smap2);
xs4_sm = sense_recon(ys4_sm, smap4);


figure(3), clf
im plc 2 3
im(1, 'row', nc, smap), cbar
im(2, 'row', nc2, smap2), cbar
im(3, 'row', nc4, smap4), cbar
% im(3, 'row', nc2, ys2(:, :, :, 1), 'ys2 fr1'), cbar
im(4, 'row', M, xs_sm), cbar
im(5, 'row', M, xs2_sm), cbar
im(6, 'row', M, xs4_sm), cbar

disp('end of test 2.  Want to see anything else before returning?')
keyboard


end


function ys = impulsedata(smap, M)

[ny nz nc] = size(smap);
xtrue = zeros(ny, nz);
impy = 3*ny/4+1;
impz = nz/2+1;
xtrue(impy, impz) = 1;

y = zeros(ny, nz, nc);
for ic = 1:nc
    sxt = smap(:,:,ic) .* xtrue;
    y(:,:,ic) = fftshift(fft2(fftshift(sxt)));
end

ys1 = zeros(size(y));
ys1(1:2:end, :, :) = y(1:2:end, :, :);
ys = repmat(ys1, [1 1 1 M]);

end
