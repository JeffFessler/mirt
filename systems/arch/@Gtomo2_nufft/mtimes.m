 function y = mtimes(ob, x)
%function y = mtimes(ob, x)	for Gtomo2_nufft object
% y = G * x, or x = G' * y

if ob.apower ~= 1, error notdone, end

%
% partial projection or backprojection (for ordered subsets)
%
if ob.is.subref
	error 'subref not done'

%
% full "projection"
%
elseif ~ob.is_transpose
	if ob.is_masked
		x = embed(x, ob.mask);
	end

	%
	% full dsft (doesn't have nxy_shift implemented, built into tomo_filter)
	%
	if ob.is.dsft
		x = reshape(x, ob.nx, ob.ny);
		y = dtft_mex('forward', ob.omega', x, int32(ob.nthread));
		y = ob.pixel_size * y;
	%	y = ob.dsft * x(:);

	%
	% nufft
	%
	else
		x = reshape(x, ob.nx, ob.ny);
		y = nufft(x, ob.st);
	end

	%
	% parallel beam case: filter, then inverse DFT
	%
	if streq(ob.geometry, 'par') % parallel
		y = spectral_filter(ob, y);

		K = ob.Krho;
		y = y([(K/2+1):K, 1:(K/2)], :);		% fft1shift
		y = ifft(y);				% [K,na] inverse 1D FFTs
		if K ~= ob.nb
			error 'over-sampling radially not done'
		end


	%----------------------------------------
	% fan beam geometry
	%----------------------------------------
	elseif streq(ob.geometry, 'fan') % fan
		y = spectral_filter(ob, y);

		% using 1Dnufft to compute 1Dnuifft
		y = y * ob.fan.del_rho;
		y = nufft(y, ob.fan.r_st); % from Krho to nb

		% interp1 along beta to get values at desired beta locations
		% Note: in fractional_delay, the delay is unitless ---the
		% fraction of sample numbers, how many samples to shift
		delay = ob.fan.sigma*(-1) / ob.fan.del_beta;
		y = fractional_delay(y.', delay.'); % delay for each column
		y = y.';

	else
		error 'Geometry not done'
	end

%
% full "back-projection"
%
else
	x = reshape(x, ob.nb, ob.na);		% sinogram shape

	K = ob.Krho;
	if K ~= ob.nb
		error 'over-sampling radially not done, need padding'
	end

	if streq(ob.geometry, 'par')
		y = (1/K) * fft(x);			% adjoint of ifft()
		y = y([(K/2+1):K, 1:(K/2)], :);		% fft1shift

	elseif streq(ob.geometry, 'fan')

		%
		% take fractional_delay for each sigma along beta, back to
		% uniform \theta (values at uniform \theta grid)
		%
		delay = ob.fan.sigma / ob.fan.del_beta;	% fix: precompute!
		y = fractional_delay(x.', delay.'); % delay for each column
		y = y.';

		%
		% adjoint of 1D NUFFT along each col
		%
		y = nufft_adj(y, ob.fan.r_st);
		y = y * conj(ob.fan.del_rho); % fix: build into tomo_filter!
	else
		error 'Geometry not done'
	end

	%
	% trick: to fix imaginary part.  subtle!
	% frankly it is a bit mysterious why this works
	% since "real()" is nonlinear so shouldn't have an adjoint!
	%
	if ob.is.shift0
		y(1,:) = 2 * real(y(1,:));
	end
	y = y .* conj(ob.tomo_filter);
	y = y(:);

	% finally, back to 2D object space
	if ob.is.dsft
		if isreal(y)
			y = complex(y);
			warning 'faking sinogram complex'
		end
		y = dtft_mex('adjoint', ob.omega', y, ...
			int32([ob.nx ob.ny]'), int32(ob.nthread));
		y = ob.pixel_size * y;
	%	y = ob.dsft' * y;

	else
		y = nufft_adj(y, ob.st);
		y = reshape(y, ob.nx, ob.ny);

	end

	if ob.is_masked
		y = y(ob.mask);
	end
end

if ~ob.is.complex
	y = real(y);
end

y = y(:);


%
% spectral_filter()
%
function y = spectral_filter(ob, y)
K = ob.Krho;
y = reshape(y, K, ob.na); % for no oversampled case in radical \rho
y = y .* ob.tomo_filter;

%
% trick: fix imaginary part.  subtle!
% see technical report for explanation of "imaginary part fix"
% x real--> y should be hermitian symmetric G(s)=conj(G(-s))
%
if ob.is.shift0
	y(1,:) = 2 * real(y(1,:));
else
	warning 'do i really not need this real trick?'
end
