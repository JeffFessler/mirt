 function y = Gtomo_nufft_filter(omega, ob, do_phase)
%function y = Gtomo_nufft_filter(omega, ob, do_phase)
%|
%| build the sinogram-spectrum-sized matrix
%| that is .* multiplied after 2D FT, before iFFT or NUiFFT
%|
%| in
%|	omega	[M 2]	frequency sample locations (radians)
%| out
%|	y	[M 2]	filter
%|
%| Copyright 2001-1, Jeff Fessler, University of Michigan
%| Extend to fan-beam geometry, 2003-11, Yingying Zhang

if nargin < 1, help(mfilename), error(mfilename), end

% effect of image-domain shift
y = exp(1i * (omega * ob.nxy_shift(:)));

% effect of image basis function (extension to non-rect by S. Matej)
% basis: b(x/Dx,y/Dy) <-FT-> Dx*Dy * B(Dx*u,Dy*v)
y = y .* (ob.dx).^2;
M = size(omega, 1);

if isempty(ob.basis.type) || streq(ob.basis.type, 'pixel')
	% sinc_2 due to square pixels
	if ob.chat, printf('pixel basis, dx=%g', ob.dx), end
	y = y .* nufft_sinc(omega(:,1)/(2*pi)) .* nufft_sinc(omega(:,2)/(2*pi));

elseif streq(ob.basis.type, 'no')
	if ob.chat, printf('no image basis modeled'), end

elseif streq(ob.basis.type, 'KB')
	if ob.chat, printf('KB basis: J=%g, alpha=%g, m=%g, n=%g', ...
		ob.basis.diam, ob.basis.shape, ob.basis.m, ob.basis.dim); end
	norm = kaiser_bessel_ft(0, ...
		ob.basis.diam, ob.basis.shape, ob.basis.m, ob.basis.dim);
	nd = sqrt(omega(:,1)/(2*pi) .* omega(:,1)/(2*pi) + ...
		omega(:,2)/(2*pi) .* omega(:,2)/(2*pi));
	y = y .* kaiser_bessel_ft(nd, ...
		ob.basis.diam, ob.basis.shape, ob.basis.m, ob.basis.dim) / norm;

elseif streq(ob.basis.type, 'Gauss')
	error('Gaussian basis not yet implemented')

else
	error(sprintf('basis function %s not implemented', ob.basis.type))
end


% form into half-sinogram shape
K = ob.Krho;
if (M == 1) && isequal(omega, [0 0]) % DC part

else
	y = reshapee(y, [], ob.na); % fix: or na/2 ?
%	if mod(K, 2) == 0
%	if size(y, 1) == K/2 * ob.na
%		y = reshape(y, [K/2 ob.na]);
%	elseif size(y, 1) == (K + 1) * ob.na
%		y = reshape(y, K + 1, ob.na);
%	end
%	else
%	if size(y, 1) == (K - 1)/2 * ob.na
%		y = reshape(y, [(K - 1)/2 ob.na]);
%	elseif size(y, 1) == K * ob.na
%		y = reshape(y, K, ob.na);
%	end
%	end
end


% detector radial response resolution-loss effect in frequency domain
% extension to parallel-beam non-rect case by S. Matej
if (M == 1) && isequal(omega, [0 0]) % DC part
	kk = 0;
else
	if mod(K, 2) == 1
		if size(y, 1) == (K - 1)/2
			kk = [1 : (K-1)/2]';
		elseif size(y, 1) == K
			kk = [- (K-1)/2 : (K-1)/2]';
		end
	else
		if size(y, 1) == K/2
			kk = [1 : K/2]';
		elseif size(y, 1) == K + 1
			kk = [-K/2 : K/2]';
		else
			kk = [-K/2 : K/2-1]'; % usual parallel-beam ?
		end
	end
end

if isempty(ob.beam.type) || streq(ob.beam.type, 'rect')

	% Trick: in fan-beam, the detector blur is not shift-invariant.
	% We approximate by the blur at the center of object.
	% The effective strip_width / ds at the center is the same.

	strip_ray = ob.strip_width / ob.ds;
	blur = nufft_sinc(strip_ray * kk / K);
	if ob.chat
		printf('strip integrals, relative beam width=%g', strip_ray)
	end

elseif streq(ob.beam.type, 'line')
	blur = ones(size(kk));
	if ob.chat, printf('line integrals modeled'), end

elseif streq(ob.beam.type, 'KB')
	if ob.chat, printf('KB beam shape: J=%g, alpha=%g, m=%g', ...
		ob.beam.diam, ob.beam.shape, ob.beam.m), end
	rel_bdiam = ob.beam.diam;	% KB_diameter relative to ob.ds
	norm = kaiser_bessel_ft(0, rel_bdiam, ob.beam.shape, ob.beam.m, 1);
	blur = kaiser_bessel_ft(...
		kk/K, rel_bdiam, ob.beam.shape, ob.beam.m, 1) / norm;

elseif streq(ob.beam.type, 'Gauss')
	g_sigma = ob.beam.shape / sqrt(8*log(2));
	blur = exp(-(2*pi*kk/K).^2 .* g_sigma^2/2);

else
	fail('beam shape %s not implemented', ob.beam.type)

end

if (M == 1) && isequal(omega, [0 0]) % DC part
	y = y .* blur;
else
	y = y .* repmat(blur, [1 ob.na]); % include blur effect
end

% phase "shift" to effect a half-pixel shift in each row
% corresponding to "offset=0" in tomographic projection.
% build in the post-fft shift too while at it.
%if streq(ob.geometry, 'par') % do we need this in fan-beam ???
if do_phase % do we need this in fan-beam ???
	phase = exp(1i*2*pi*(K+ob.is.shift0)/2 * kk / K);
	y = y .* repmat(phase, [1 ob.na]);
	y = y ./ ob.ds;	% see JF tech. report
end

% trick: for parallel case, build in the phase shift due to offset_s
% added 2005-8-24 since it had been omitted previously
if ~ob.is.fan
	phase = exp(-2i*pi * ob.offset_s * kk / K);
	y = y .* repmat(phase, [1 ob.na]);
end
