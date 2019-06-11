 function ob = Gtomo_blob(mask, dim_data, varargin)
%function ob = Gtomo_blob(mask, dim_data, options)
%|
%| Construct Gtomo_blob object, for fast forward and backprojection.
%| INCOMPLETE, WORK-IN-PROGRESS!
%|
%| Currently handles 2D parallel and fan-beam cases.
%|
%| See Gtomo_blob_test.m for example usage.
%|
%| in
%|	mask	[nx,ny]	logical array of object support.
%|	dim_data [2]	data dimensions, usually: [nb na], where
%|		nb	number of "detector bins" in each sinogram row
%|		na	number of view angles: size(sinogram,2)
%|
%| options
%|	'chat'		verbose printing of debug messages
%|	'pixel_size'	width of pixels
%|	'ray_spacing'	radial sample spacing
%|	'strip_width'	width of rectangular detector PSF
%|	'orbit'		projection angle coverage [180 degrees]
%|	'orbit_start'	first projection angle [0 degrees]
%|	'xscale'	use -1 to flip in x direction
%|	'yscale'	use -1 to flip in y direction
%|
%|	'beam'		structure for beam PSF (see below)
%|	'basis'		structure for image-domain basis (see below)
%|
%| required arguments for fan-beam geometry (for which a value must follow):
%|	'dis_src_det'		source-to-detector distance
%|	'dis_iso_det'		object_isocenter-to-detector distance
%|	'source_offset'		offset between source-detector and isocenter
%|	'channel_offset'	offset between central ray and detector center
%|				in terms of fraction of sample(detector) spacing
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|
%| For more help, see Gtomo_blob_test.m
%|
%| Copyright 2005-7-26, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mask, 'test'), Gtomo_blob_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;
[arg.nx arg.ny] = size(mask);
arg.nb = dim_data(1);
arg.na = dim_data(2);

% option defaults
arg.chat = 1;
arg.nthread = 1;

arg.pixel_size = 1;
arg.ray_spacing = 1;
arg.strip_width = 1;
arg.orbit = 360;
arg.orbit_start = 0;
arg.xscale = 1;
arg.yscale = 1;

arg.dis_src_det = [];
arg.dis_iso_det = [];
arg.source_offset = 0;
arg.channel_offset = 0;

%
% image basis function used in the iterative reconstruction model (S. Matej)
% types:
%	'pixel'
%	'KB'	blob: diam=J, shape=alpha, m=KB_order, dim=dimension (2D/3D),
%		kernel - blob values (image=kernel*img_coefs, * is convolution)
%	'Gauss'	Gaussian - shape=FWHM (diam=diameter of truncation)
%	'no'	no basis function modeled/considered in the blob
%
arg.basis.type	= 'kb';
arg.basis.dim	= ndims(mask);
arg.basis.diam	= 4;	% diameter relative to pixel_size
arg.basis.m	= 2;
arg.basis.kernel = [];	% JxJ matrix
arg.basis.shape	= 10.4;	% 2d default

% model for beam radial PSF
arg.beam.type	= 'rect';	% 'rect', 'line', 'Gauss', 'KB'
arg.beam.diam	= 1;		% width of Gauss and KB relative to ray_spacing
arg.beam.shape	= 0;
arg.beam.m	= 0;

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

% initialize
arg = Gtomo_blob_setup(arg);

arg.nd = prod(dim_data);
arg.np = sum(mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtomo_blob', ...
	'forw', @Gtomo_blob_forw, 'back', @Gtomo_blob_back);


%
% Gtomo_blob_setup()
%
function arg = Gtomo_blob_setup(arg)

arg.is.fan = ~isempty(arg.dis_src_det);
if arg.is.fan % trick: isinf([]) returns []
	arg.is.fan = arg.is.fan & ~isinf(arg.dis_src_det);
end

if arg.is.fan
	arg = Gtomo_blob_init_fan(arg);
else
	arg = Gtomo_blob_init_par(arg);
end

% fix: need "center_x and center_y" here, and consider flip_y or yscale

%if arg.chat, disp 'init_basis', end
%arg = Gtomo_blob_init_basis(arg);	% space-domain basis function


%
% Gtomo_blob_init_par()
%
function arg = Gtomo_blob_init_par(arg)

k_ang = ([0:(arg.na-1)]/arg.na * arg.orbit + arg.orbit_start) * pi / 180;
k_rad = [-arg.Krho/2:arg.Krho/2-1] * (2*pi/arg.Krho) * arg.pixel_size / arg.ray_spacing;
[rr aa] = ndgrid(k_rad, k_ang);

% desired k-space coordinates in [-pi,pi)
arg.omega = [arg.xscale * col(rr .* cos(aa)), -arg.yscale * col(rr .* sin(aa))];

% fix: should check withing mask only!
if max(arg.nx,arg.ny) * arg.pixel_size > arg.nb * arg.ray_spacing
	error 'truncated projections probably will not work now'
end
if arg.pixel_size > arg.ray_spacing
	warning 'big pixels and small rays may not work without oversampling k'
end


%
% Gtomo_blob_init_fan()
%
function arg = Gtomo_blob_init_fan(arg)

arg.dis_src_iso = arg.dis_src_det - arg.dis_iso_det;

% check to see if object lies within scanner FOV
fov_fan_check(arg)

%-------------------------------------------------------
% finding the desired locations in fan(na, nb) data
%-------------------------------------------------------
fan.beta = ([0:(arg.na-1)]/arg.na * arg.orbit + arg.orbit_start) * pi / 180;
fan.del_beta = diff(fan.beta([1 2]));

fan.del_sigma = arg.ray_spacing / arg.dis_src_det; % in radians
fan.sigma = ([-(arg.nb-1)/2:(arg.nb-1)/2]' ...
			- arg.channel_offset) * fan.del_sigma;

% Note: in fractional_delay, the delay is unitless sample fractions.
arg.delay1 = fan.sigma' / fan.del_beta;	% precompute angular shift

%-------------------------------------------------------------------
% find desired location in R (corresponding to sigma) for 1D NUIFFT
% (@each angle: non in rad)
%-------------------------------------------------------------------

% decide from relation in PENG paper
fan.dis_src_iso = arg.dis_src_det - arg.dis_iso_det;
fan.sin_rad = fan.dis_src_iso * sin(fan.sigma);

if isempty(arg.rmax)
	% fix: this next line looks wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	arg.rmax = abs(fan.sin_rad(end) - fan.sin_rad(1));
end
arg.fan.del_rho = 1/arg.rmax; % it has to be del_rho <= 1/rmax
arg.fan.omegam = (-1) * 2*pi * arg.fan.del_rho * fan.sin_rad;

% center shift shouldn't be 0 since it wasn't absorbed into tomo_filter
int1c = arg.interp; % copy for 1d
if ~isempty(arg.interp) && streq(arg.interp{1}, 'table')
	interp1 = {arg.interp{end}};
	if isnumeric(interp1{:})
		int1c = {'minmax:kb'};
%		arg.fan.r_st = blob_init(arg.fan.omegam, arg.Krho, arg.J, arg.Krho*2, arg.Krho/2, 'minmax:kb');
	elseif ischar(interp1{:})
		int1c = interp1;
%		arg.fan.r_st = blob_init(arg.fan.omegam, arg.Krho, arg.J, arg.Krho*2, arg.Krho/2, interp1{:});
	else
		error 'wrong interp input'
	end
%else
%	arg.fan.r_st = blob_init(arg.fan.omegam, arg.Krho, arg.J, arg.Krho*2, arg.Krho/2, arg.interp{:});
end

if arg.is.dsft1
	% no need for 1D blob structure
	arg.fan.n_shift = arg.Krho / 2;
else
	arg.fan.r_st = blob_init(arg.fan.omegam, ...
		arg.Krho, arg.J, arg.Krho*2, arg.Krho/2, int1c{:});
end

%------------------------------------------------------------------
% find corresponding k-space polar locations (uniform in both)
%------------------------------------------------------------------
k_ang = fan.beta; % \theta when sigma=0

%del_k_rad = 2*pi / arg.Krho;
k_rad = 2*pi * [-arg.Krho/2:arg.Krho/2-1]' * arg.fan.del_rho * arg.pixel_size;
if arg.chat
	printf('k_rad size %g %g k_ang size %g %g', ...
		size(k_rad), size(k_ang))
end
[rr aa] = ndgrid(k_rad, k_ang);

% desired k-space coordinates in [-pi,pi)
arg.omega = [arg.xscale * col(rr .* cos(aa)), -arg.yscale * col(rr .* sin(aa))];

if arg.pixel_size > arg.ray_spacing
	warning 'big pixels and small rays may require oversampling k'
end


%
% fov_fan_check()
%
function fov_fan_check(arg)
x = ([0:arg.nx-1] - (arg.nx-1)/2) * arg.pixel_size;
y = ([0:arg.ny-1] - (arg.ny-1)/2) * arg.pixel_size;
[x y] = ndgrid(x, y);
r = sqrt(x.^2 + y.^2);
%fov_max = sqrt(nx^2 + ny^2) / 2 * arg.pixel_size;
mask_max = max(r(arg.mask(:)));
s_max = ((arg.nb-1)/2 - abs(arg.channel_offset)) * arg.ray_spacing;
r_max = arg.dis_src_iso * sin(s_max / arg.dis_src_det); % scanner FOV radius
if mask_max > r_max
	printf('Scanner FOV radius = %g', r_max)
	printf('Object support maximum radius = %g', mask_max)
	printf 'Object support exceeds FOV of this system'
	warning 'Truncated projections probably will not work, as implemented'
end


%
% Gtomo_blob_init_basis()
% initialization of the image basis function (blob/Gauss/...)
% by S. Matej
%
function ob = Gtomo_blob_init_basis(ob)

if streq(ob.basis.type, 'KB')
	dim1 = floor(ob.basis.diam/2);
	dim2 = 1+ 2*dim1;
	row = [-dim1:dim1] .* [-dim1:dim1];
	onev(1:dim2) = 1;
	kernel = row'*onev + onev'*row;
	%normft = ...
	%kaiser_bessel_ft(0,ob.basis.diam ,ob.basis.shape,ob.basis.m,ob.basis.dim)
	if ob.basis.dim==2
		kernel = sqrt(kernel);
		ob.basis.kernel = kaiser_bessel( ...
			kernel, ob.basis.diam, ob.basis.shape, ob.basis.m);

	elseif ob.basis.dim==3
		for i = 1:dim2
			kernel3(:,:,i) = kernel+row(i);
		end
		kernel3 = sqrt(kernel3);
		ob.basis.kernel = kaiser_bessel(...
			kernel3, ob.basis.diam, ob.basis.shape, ob.basis.m);

	else
		fail('basis function dimension %g not supported', ob.basis.dim)
	end

	norm=sum(sum(sum(ob.basis.kernel)));
	ob.basis.kernel = ob.basis.kernel/norm;

elseif streq(ob.basis.type,'Gauss')
	error('Gaussian basis not yet implemented - kernel not calculated')

elseif isempty(ob.basis.type) || streq(ob.basis.type,'pixel') || ...
	streq(ob.basis.type,'no')

else
	fail('basis function %s not implemented', ob.basis.type)

end


%
% Gtomo_blob_forw(): y = G * x
%
function y = Gtomo_blob_forw(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np);

%
% full dsft (doesn't have nxy_shift implemented, built into tomo_filter)
%
if arg.is.dsft2
	y = jf_mex('dtft,forward', arg.omega', x, int32(arg.nthread));

%
% blob
%
else
	y = blob(x, arg.st);
end

y = Gtomo_blob_spectral_filter(arg, y);

%
% parallel beam case: filter, then inverse DFT
%
if ~arg.is.fan % parallel

	K = arg.Krho;
%	y = y([(K/2+1):K, 1:(K/2)], :);		% fft1shift
	y = fftshift(y, 1);		% fft1shift
	y = ifft(y);				% [K,na] inverse 1D FFTs
	if K > arg.nb
		y = y(1:arg.nb,:);
		warning 'over-sampling radially not tested'
	end


%
% fan beam
%
else
	% using 1D blob to compute 1D nuifft
	if arg.is.dsft1
		y = jf_mex('dtft,forward', arg.fan.omegam', y, int32(arg.nthread));
		% phase effect due to n_shift:
		y = y .* repmat(exp(1i*arg.fan.omegam*arg.fan.n_shift), [1 ncol(y)]);
	else
		y = blob(y, arg.fan.r_st); % from Krho to nb
	end

	% interp1 along beta to get values at desired beta locations
	y = fractional_delay(y.', -arg.delay1).'; % delay for each column
end

if ~arg.is.complex
	y = real(y);
end

y = ei.shape(y);


%
% Gtomo_blob_back(): x = G' * y
% full backprojection
%
function x = Gtomo_blob_back(arg, y)

[y eo] = embed_out(y, [arg.nb arg.na]);

K = arg.Krho;
if K > arg.nb
	y = [y; zeros(K-arg.nb, arg.na)];
	warning 'over-sampling radially not tested and probably wrong!'
end

if ~arg.is.fan
	y = (1/K) * fft(y);			% adjoint of ifft()
%	x = x([(K/2+1):K, 1:(K/2)], :);		% fft1shift
	y = fftshift(y, 1);	% fft1shift

else

	%
	% take fractional_delay for each sigma along beta, back to
	% uniform \theta (values at uniform \theta grid)
	%
	y = fractional_delay(y.', arg.delay1).'; % delay for each column

	%
	% adjoint of 1D NUFFT along each col
	%
	if arg.is.dsft1
		% phase effect due to n_shift:
		y = y .* repmat(exp(-1i*arg.fan.omegam*arg.fan.n_shift), [1 ncol(y)]);
		y = jf_mex('dtft,adjoint', arg.fan.omegam', y, ...
			int32(arg.Krho), int32(arg.nthread));
	else
		y = blob_adj(y, arg.fan.r_st);
	end
end

%
% trick: to fix imaginary part.  subtle!
% frankly it is a bit mysterious why this works
% since "real()" is nonlinear so shouldn't have an adjoint!
%
if arg.is.shift0
	y(1,:) = 2 * real(y(1,:));
end
y = y .* conj(arg.tomo_filter);
y = y(:);

% finally, back to 2D object space
if arg.is.dsft2
	if isreal(y)
		y = complex(y);
		warning 'faking sinogram complex'
	end
	x = jf_mex('dtft,adjoint', arg.omega', y, ...
			int32([arg.nx arg.ny]'), int32(arg.nthread));

else
	x = blob_adj(y, arg.st);
	x = reshape(x, arg.nx, arg.ny);
end

if ~arg.is.complex
	x = real(x);
end

x = eo.shape(x, arg.mask, arg.np);


%
% Gtomo_blob_spectral_filter()
%
function y = Gtomo_blob_spectral_filter(ob, y)
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
