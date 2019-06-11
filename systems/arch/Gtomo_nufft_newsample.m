 function ob = Gtomo_nufft_newsample(mask, dim_data, varargin)
%function ob = Gtomo_nufft_newsample(mask, dim_data, options)
%
% Construct Gtomo_nufft object, for Fourier-based forward and backprojection.
% Currently handles 2D parallel and fan-beam cases.
%
% See Gtomo_nufft_test.m for example usage.
%
% Basically, you create an object calling:
%		G = Gtomo_nufft(...)
% and then you can use it thereafter by typing commands like
%		y = G * x;
% which will auto-magically do the multiplication by calling a mex file.
%
% in
%	mask	[nx,ny]	logical array of object support.
%	dim_data [2]	data dimensions, usually: [nb na], where
%		nb	number of "detector bins" in each sinogram row
%		na	number of view angles: size(sinogram,2)
%
% optional arguments for which a value must follow:
%	'chat'		verbose printing of debug messages
%	'pixel_size'	width of pixels
%	'ray_spacing'	radial sample spacing
%	'strip_width'	width of rectangular detector PSF
%	'orbit'		projection angle coverage [180 degrees]
%	'orbit_start'	first projection angle [0 degrees]
%	'xscale'	use -1 to flip in x direction
%	'yscale'	use -1 to flip in y direction
%
%	'beam'		structure for beam PSF (see below)
%	'basis'		structure for image-domain basis (see below)
%
%	'is.dsft'	use full DSFT rather than NUFFT
%			(this is expensive but useful for debugging)
%	'is.shift0'	apply shifts so center of rotation is at 0 (is default)
%	'is.complex'	return complex projection and adjoint (unusual)
%
% required arguments for fan-beam geometry (for which a value must follow):
%	'dis_src_det'		source-to-detector distance
%	'dis_iso_det'		object_isocenter-to-detector distance
%	'source_offset'		offset between source-detector and isocenter
%	'channel_offset'	offset between central ray and detector center
%				in terms of fraction of sample(detector) spacing
%
% interpolation specification
%
%	'J'		neighborhood size (default is 6)
%	'Kd'		# of FFT points along each dimension 2*[nx ny]
%	'Krho'		number of radial samples in polar space, default is nb.
%
%	The default (recommended) interpolator uses minmax with KB scaling.
%	These others below are provided for making comparisons.
%	'is.kaiser'	use minmax tuned kaiser-bessel interpolator
%	'interp'	{ cell array of arguments to pass to nufft_init() }
%			(This allows complete user control of interpolator.)
%	Example: straw-man bilinear interpolation is:
%		..., 'interp', {'linear'}, 'J', 2, ...
%	Example: for straw-man uniform scaling factors in NUFFT:
%		..., 'interp', {'minmax:user', {{1}, {1}}, {{0.5}, {0.5}}}
%	Example: for table-based NUFFT for large-scale problems:
%		..., 'interp', {'table', 'minmax:kb'}
%
% out
%	ob [nd,np]	np = sum(mask(:)), so it is already "masked"
%
% For more help, see Gtomo_nufft_test.m
%
% Copyright 2005-1-19, Jeff Fessler, The University of Michigan
% Fan-beam contributed by Yingying Zhang.  Basis etc. by Samuel Matej.

if nargin == 1 && streq(mask, 'test'), Gtomo_nufft_test, return, end
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

% fix: we could over-sample in the future, especially if nb not a power of 2!
%arg.Krho = 2*ceil(arg.nb/2);
arg.rmax = []; % maximum object radius in FOV, relevant to sampling
arg.Krho = arg.nb;
arg.interp = {};
arg.is.kaiser = false;
arg.J = 6;
arg.Kd = 2 * size(mask);

arg.is.dsft = false;
arg.is.complex = false;
arg.is.shift0 = true;

%
% image basis function used in the iterative reconstruction model (S. Matej)
% types:
%	'pixel'
%	'KB'	blob: diam=J, shape=alpha, m=KB_order, dim=dimension (2D/3D),
%		kernel - blob values (image=kernel*img_coefs, * is convolution)
%	'Gauss'	Gaussian - shape=FWHM (diam=diameter of truncation)
%	'no'	no basis function modeled/considered in the nufft
%
arg.basis.type	= 'pixel';
arg.basis.diam	= 1;	% Relative_diameter - relative to pixel_size
arg.basis.shape	= 0;
arg.basis.m	= 0;
arg.basis.dim	= 2;
arg.basis.kernel = [];	% JxJ matrix

% model for beam radial PSF
arg.beam.type	= 'rect';	% 'rect', 'line', 'Gauss', 'KB'
arg.beam.diam	= 1;		% width of Gauss and KB relative to ray_spacing
arg.beam.shape	= 0;
arg.beam.m	= 0;

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

% initialize geometry stuff, parallel or fan
arg = Gtomo_nufft_setup(arg);

arg.nd = prod(dim_data);
arg.np = sum(mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtomo_nufft', ...
	'forw', @Gtomo_nufft_forw, 'back', @Gtomo_nufft_back);


%
% Gtomo_nufft_setup()
%
function arg = Gtomo_nufft_setup(arg)

arg.is.fan = ~isempty(arg.dis_src_det);
if arg.is.fan % trick: isinf([]) returns []
	arg.isfan = arg.is.fan & ~isinf(arg.dis_src_det);
end

if arg.is.fan
	arg = Gtomo_nufft_init_fan(arg);
else
	arg = Gtomo_nufft_init_par(arg);
end

% fix: need "center_x and center_y" here, and consider flip_y or yscale
arg.nxy_shift = size(arg.mask)/2 - arg.is.shift0 * 0.5;	% shifts

%
% create the structure needed by NUFFT
%
if arg.is.dsft
	% no need to precompute anything
% arg.dsft = dsft2_init(arg.omega, nx, ny, 0 !! & arg.is.shift0, arg.yscale);

else % precompute for NUFFT

	if arg.is.kaiser
		arg.interp = {'kaiser'};
	elseif isempty(arg.interp)
		arg.interp = {'minmax:kb'};	% minmax based on KB scaling
	end

	% trick: made shift 0 since built into tomo_filter!
	if arg.chat, disp 'nufft_init', end
	arg.st = nufft_init(arg.omega, [arg.nx arg.ny], [arg.J arg.J], arg.Kd, ...
		0 * arg.nxy_shift, arg.interp{:});
end

if arg.chat, disp 'init_basis', end
arg = Gtomo_nufft_init_basis(arg);	% space-domain basis function

%
% make filter to be used in sinogram spectral domain
%
if ~arg.is.fan
    arg.tomo_filter = Gtomo_nufft_filter(arg.omega, arg, ~arg.is.fan);
else
    % arg.omega is for (\rho > 0): no DC
    arg.tomo_filter = Gtomo_nufft_filter_yy(arg.omega, arg, ~arg.is.fan);
	arg.tomo_filter = arg.tomo_filter * arg.fan.del_rho;
end


%
% Gtomo_nufft_init_par()
%
function arg = Gtomo_nufft_init_par(arg)

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
% Gtomo_nufft_init_fan()
%
function arg = Gtomo_nufft_init_fan(arg)

arg.dis_src_iso = arg.dis_src_det - arg.dis_iso_det;

% check to see if object lies within scanner FOV
fov_fan_check(arg)

if arg.orbit ~= 360
	error 'only 360 orbit implemented due to periodic fractional delay'
end

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
    % the way calculating del_rho makes k_rad is outside of [-pi pi]
	arg.rmax = abs(fan.sin_rad(end) - fan.sin_rad(1)) + abs(fan.sin_rad(2) - fan.sin_rad(1)) ;
end
arg.fan.del_rho = 1/arg.rmax; % it has to be del_rho <= 1/rmax
arg.fan.omegam = (-1) * 2*pi * arg.fan.del_rho * fan.sin_rad;

% center shift shouldn't be 0 since it wasn't absorbed into tomo_filter
int1c = arg.interp;
if ~isempty(arg.interp) & streq(arg.interp{1}, 'table')
	interp1 = {arg.interp{end}};
	if isnumeric(interp1{:})
		int1c = {'minmax:kb'};
	elseif ischar(interp1{:})
		int1c = interp1;
	else
		error 'wrong interp input'
	end
end

%------------------------------------------------------------------
% find corresponding k-space polar locations (uniform in both)
%------------------------------------------------------------------
k_ang = fan.beta; % \theta when sigma=0

if mod(arg.Krho, 2) == 1
    k_rad = 2*pi * [-(arg.Krho - 1)/2:(arg.Krho - 1)/2]' * arg.fan.del_rho * arg.pixel_size;
elseif mod(arg.Krho, 2) == 0
    k_rad = 2*pi * [-arg.Krho/2:arg.Krho/2]' * arg.fan.del_rho * arg.pixel_size;
end

Kr = size(k_rad, 1);

if arg.chat
	printf('k_rad size %g %g k_ang size %g %g', ...
		size(k_rad), size(k_ang))
end
[rr1 aa1] = ndgrid(k_rad, k_ang);

[rr aa] = ndgrid(k_rad((Kr - 1)/2 + 2:end), k_ang);

% desired k-space coordinates in [-pi,pi)
% arg.allomega = [arg.xscale * col(rr1 .* cos(aa1)), -arg.yscale * col(rr1 .* sin(aa1))];
arg.omega = [arg.xscale * col(rr .* cos(aa)), -arg.yscale * col(rr .* sin(aa))];

arg.fan.r_st = ...
	nufft_init(arg.fan.omegam, (Kr-1)/2, 6, Kr-1, -1, int1c{:}); % shift acct for index from 0, not 1

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
% Gtomo_nufft_init_basis()
% initialization of the image basis function (blob/Gauss/...)
% by S. Matej
%
function ob = Gtomo_nufft_init_basis(ob)

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
		error(sprintf('basis function dimension %g not supported', ob.basis.dim))
	end

	norm=sum(sum(sum(ob.basis.kernel)));
	ob.basis.kernel = ob.basis.kernel/norm;

elseif streq(ob.basis.type,'Gauss')
	error('Gaussian basis not yet implemented - kernel not calculated')

elseif isempty(ob.basis.type) | streq(ob.basis.type,'pixel') | ...
	streq(ob.basis.type,'no')

else
	error(sprintf('basis function %s not implemented', ob.basis.type))

end


%
% Gtomo_nufft_forw(): y = G * x
%
function y = Gtomo_nufft_forw(arg, x)

% if needed, convert concise column to 3d array
flag_column = 0;
if size(x,1) == arg.np
	flag_column = 1;
	x = embed(x, arg.mask);
end

%
% full dsft (doesn't have nxy_shift implemented, built into tomo_filter)
%
if arg.is.dsft
	y = dtft_mex('forward', arg.omega', x, int32(arg.nthread));
	y = arg.pixel_size * y;
%	y = arg.dsft * x(:);

%
% nufft
%
else
    y = nufft(x, arg.st);
end

% for DC
y0 = sum(x(:)); % dc
y0 = y0 .* Gtomo_nufft_filter_yy([0 0], arg, ~arg.is.fan) .* arg.fan.del_rho;

y = Gtomo_nufft_spectral_filter(arg, y);

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
	% using 1D nufft to compute 1D nuifft

	y = nufft(y, arg.fan.r_st); % from Krho/Krho+1 to nb
	y = 2 .* real(y) + ones(size(y)) .* y0;

	% interp1 along beta to get values at desired beta locations
	y = fractional_delay(y.', -arg.delay1).'; % delay for each column

    % use linear interp to see if there is circular err pattern
%     for ib = 1:arg.nb
%         y(ib,:) = perodic_shift(y(ib,:), arg.delay1(ib));
%     end
end

if ~arg.is.complex
	y = real(y);
end

if flag_column % column in yields column out.
	y = y(:);
end


%
% Gtomo_nufft_back(): x = G' * y
% full backprojection
%
function x = Gtomo_nufft_back(arg, y)

flag_column = 0;
if size(y,1) == arg.nd
	flag_column = 1;
	y = reshape(y, [arg.nb arg.na]); % fix: what if multiple columns?
end

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
    y0 = sum(y(:)) .* conj(Gtomo_nufft_filter_yy([0 0], arg, ~arg.is.fan) .* arg.fan.del_rho);

%      for ib = 1:arg.nb
%         y(ib,:) = perodic_shift(y(ib,:), -arg.delay1(ib));
%      end

	%
	% adjoint of 1D NUFFT along each col
	%
	y = nufft_adj(y, arg.fan.r_st);
%     y1 = [conj(flipud(y)); y0; y];
end

%
% trick: to fix imaginary part.  subtle!
% frankly it is a bit mysterious why this works
% since "real()" is nonlinear so shouldn't have an adjoint!
%
if arg.is.shift0 & ~arg.is.fan
	y(1,:) = 2 * real(y(1,:));
end
y = y .* conj(arg.tomo_filter);

% if mod(arg.Krho,2) == 0
%     y = y1(arg.Krho/2 + 2:end, :);
%     y0 = sum(y1(arg.Krho/2 + 1,:)); % dc
% else
%     y = y1((arg.Krho - 1)/2 + 2:end, :);
%     y0 = sum(y1((arg.Krho - 1)/2 + 1, :));
% end

% finally, back to 2D object space
if arg.is.dsft
	if isreal(y)
		y = complex(y);
		warning 'faking sinogram complex'
	end
	x = dtft_mex('adjoint', arg.omega', y, ...
			int32([arg.nx arg.ny]'), int32(arg.nthread));
	x = arg.pixel_size * x;
%	x = ob.dsft' * y;

else
    y = y(:);
    x = nufft_adj(y, arg.st);

%     y0 = sum(y);
%     xx = nufft_adj(y, arg.st);
    x = y0 + 2 * real(x);
	x = reshape(x, arg.nx, arg.ny);
end

if ~arg.is.complex
	x = real(x);
end

if flag_column
	x = x(arg.mask);
end


%
% Gtomo_nufft_spectral_filter()
%
function y = Gtomo_nufft_spectral_filter(ob, y)

K = ob.Krho;
% y = reshape(y, K, ob.na); % for no oversampled case in radical \rho
if mod(K, 2) == 1
    y = reshape(y, [(K - 1)/2 ob.na]);
%     y = reshape(y, K, ob.na);
else
    y = reshape(y, [K/2 ob.na]);
%     y = reshape(y, K + 1, ob.na);
end

y = y .* ob.tomo_filter;

%
% trick: fix imaginary part.  subtle!
% see technical report for explanation of "imaginary part fix"
% x real--> y should be hermitian symmetric G(s)=conj(G(-s))
%
if ob.is.shift0 & ~ob.is.fan
	y(1,:) = 2 * real(y(1,:));
else
	warning 'do i really not need this real trick?'
end
