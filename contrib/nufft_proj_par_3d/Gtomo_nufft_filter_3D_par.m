 function y = Gtomo_nufft_filter_3D_par(omega, ob, do_phase)
%function y = Gtomo_nufft_filter_3D_par(omega, ob, do_phase)
% build the sinogram-spectrum-sized matrix
% that is .* multiplied after 3D FT, before iFFT or NUiFFT
% in
%	omega	[M 3]	frequency sample locations (radians)
% out
%	y	[M 3]	filter
% Based on Gtomo_nufft_filter 
% (Copyright 2001-1, Jeff Fessler, The University of Michigan)
% (Extend to fan-beam geometry, 2003-11, Yingying Zhang)
% Extend to 3D parallel-beam geometry, 2007-04-16 Yong Long


% effect of image-domain shift
y = exp(1i * (omega * ob.nxyz_shift(:)));

%
% effect of image basis function (extension to non-rect by S. Matej)
% basis: b(x/Dx,y/Dy) <-FT-> Dx*Dy * B(Dx*u,Dy*v)
%
y = y .* (ob.dx).^3;
M = size(omega, 1);

if streq(ob.basis.type, 'cube')
	% sinc_3 due to cubes
	if ob.chat, printf('cubic basis, dx=%g', ob.dx), end
	y = y .* nufft_sinc(omega(:,1)/(2*pi)) .* nufft_sinc(omega(:,2)/(2*pi)) .* ...
        nufft_sinc(omega(:,3)/(2*pi));
else
	error(sprintf('basis function %s not implemented', ob.basis.type))
end

% form into half-sinogram shape
Kv1 = ob.Kv1;
Kv2 = ob.Kv2;
if (M == 1) & (omega == [0 0 0]) % DC part
else
    y = reshape(y, [Kv1 ob.na (Kv2+1)/2 ob.npo]);   
    % arrange data to format: (v1, v2, phi, theta)
    y = permute(y, [1 3 2 4]); 
end  

%
% detector radial response resolution-loss effect in frequency domain
%
if (M == 1) & (omega == [0 0 0]) % DC part
    kkv1 = 0;
    kkv2 = 0;
else
    if mod(Kv1, 2) == 1
        if size(y, 1) == (Kv1 - 1)/2
            kkv1 = [1 : (Kv1-1)/2]';
        elseif size(y, 1) == Kv1
            kkv1 = [- (Kv1-1)/2 : (Kv1-1)/2]';
        elseif size(y, 1) == Kv1-1
            kkv1 = [-(Kv1-1)/2 : (Kv1-1)/2-1]';
        end
    else
        if size(y, 1) == Kv1/2
            kkv1 = [1 : Kv1/2]';
        elseif size(y, 1) == Kv1 + 1
            kkv1 = [-Kv1/2 : Kv1/2]';
        else
            kkv1 = [-Kv1/2 : Kv1/2-1]'; % usual parallel-beam ?
        end
    end
    
    if mod(Kv2, 2) == 1
        if size(y, 2) == (Kv2 - 1)/2
            kkv2 = [1 : (Kv2-1)/2]';
        elseif size(y, 2) == (Kv2 + 1)/2
            kkv2 = [0 : (Kv2-1)/2]';
        elseif size(y, 2) == Kv2
            kkv2 = [- (Kv2-1)/2 : (Kv2-1)/2]';
        elseif size(y, 2) == Kv2-1
            kkv2 = [-(Kv2-1)/2 : (Kv2-1)/2-1]'; 
        end
    else
        if size(y, 2) == Kv2/2
            kkv2 = [1 : Kv2/2]';
        elseif size(y, 2) == Kv2 + 1
            kkv2 = [-Kv2/2 : Kv2/2]';
        else
            kkv2 = [-Kv2/2 : Kv2/2-1]'; % usual parallel-beam ?
        end
    end
end

%detector blur: rect2
if isempty(ob.beam.type) | streq(ob.beam.type, 'rect2')
    strip_ray_s = ob.rect_s / ob.ds;
    strip_ray_t = ob.rect_t / ob.dt;
    
    [bv1 bv2] = ndgrid(strip_ray_s * kkv1 / Kv1, strip_ray_t * kkv2 / Kv2);
    blur = nufft_sinc(bv1) .* nufft_sinc(bv2);
    
    if ob.chat
        printf('cuboid integrals, relative beam length=%g, width=%g', strip_ray_s, strip_ray_t)
    end
else 
    error(sprintf('beam shape %s not implemented', ob.beam.type))
end

if (M == 1) & (omega == [0 0 0]) % DC part
    y = y .* blur;
else
    y = y .* repmat(blur, [1  1  ob.na ob.npo]);	% include blur effect
end

%
% phase "shift" to effect a half-pixel shift in each row
% corresponding to "offset=0" in tomographic projection.
% build in the post-fft shift too while at it.
%
if do_phase 
    [kkkv1 kkkv2] = ndgrid(kkv1, kkv2);
    phase = exp(1i * 2 * pi * ((Kv1 + ob.is.shift0)/2 * kkkv1 / Kv1 ...
        + (Kv2 + ob.is.shift0)/2 * kkkv2 / Kv2));
    y = y .* repmat(phase, [1 1 ob.na ob.npo]);
    y = y ./ ob.ds ./ ob.dt;
end

% trick: for parallel case, build in the phase shift due to offset_s
% added 2005-8-24 since it had been omitted previously
if do_phase
    phase = exp(-2i*pi * (ob.offset_s * kkkv1/ Kv1 + ob.offset_t * kkkv2 / Kv2));
    y = y .* repmat(phase, [1 1 ob.na ob.npo]);
end
