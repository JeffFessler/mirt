 function ob = Gtomo_nufft(mask, dim_data, varargin)
%function ob = Gtomo_nufft(mask, dim_data, options)
%|
%| Construct Gtomo_nufft object, for Fourier-based forward and backprojection.
%| Currently handles 2D parallel and fan-beam cases.
%|
%| See Gtomo_nufft_test.m for example usage.
%|
%| Basically, you create an object calling:
%|		G = Gtomo_nufft(...)
%| and then you can use it thereafter by typing commands like
%|		y = G * x;
%| which will auto-magically do the multiplication by calling a mex file.
%|
%| in
%|	mask	[nx,ny]	logical array of object support.
%|	dim_data [2]	data dimensions, usually: [nb na], where
%|		nb	number of "detector bins" in each sinogram row
%|		na	number of view angles: size(sinogram,2)
%|
%| optional arguments for which a value must follow:
%|	'chat'			verbose printing of debug messages
%|	'dx' | 'pixel_size'	pixel width
%|	'ds' | 'ray_spacing'	radial sample spacing
%|	'strip_width'		width of rectangular detector PSF
%|	'orbit'		projection angle coverage [180 degrees]
%|	'orbit_start'	first projection angle [0 degrees]
%|	'xscale'	use -1 to flip in x direction
%|	'yscale'	use -1 to flip in y direction
%|
%|	'beam'		structure for beam PSF (see below)
%|	'basis'		structure for image-domain basis (see below)
%|
%|	'is.dsft1'	use full 1D DSFT rather than 1D NUFFT
%|	'is.dsft2'	use full 2D DSFT rather than 2D NUFFT
%|			(these are expensive but useful for debugging)
%|	'is.shift0'	apply shifts so center of rotation is at 0 (is default)
%|	'is.complex'	return complex projection and adjoint (unusual)
%|	'class'		'Fatrix' or 'fatrix2' (default)
%|
%| required arguments for fan-beam geometry (for which a value must follow):
%|	'dis_src_det'		source-to-detector distance
%|	'dis_iso_det'		object_isocenter-to-detector distance
%|	'source_offset'		offset between source-detector and isocenter
%|	'offset_s' | 'channel_offset'	offset between central ray and detector center
%|				in terms of fraction of sample(detector) spacing
%|
%| interpolation specification
%|
%|	'J'		neighborhood size (default is 6)
%|	'Kd'		# of FFT points along each dimension 2*[nx ny]
%|	'Krho'		number of radial samples in polar space, default is nb.
%|
%|	The default (recommended) interpolator uses minmax with KB scaling.
%|	These others below are provided for making comparisons.
%|	'is.kaiser'	use minmax tuned kaiser-bessel interpolator
%|	'interp'	{ cell array of arguments to pass to nufft_init() }
%|			(This allows complete user control of interpolator.)
%|	Example: straw-man bilinear interpolation is:
%|		..., 'interp', {'linear'}, 'J', 2, ...
%|	Example: for straw-man uniform scaling factors in NUFFT:
%|		..., 'interp', {'minmax:user', {{1}, {1}}, {{0.5}, {0.5}}}
%|	Example: for table-based NUFFT for large-scale problems:
%|		..., 'interp', {'table', 2^11, 'minmax:kb'}
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|
%| For more help, see Gtomo_nufft_test.m
%|
%| Copyright 2005-1-19, Jeff Fessler, University of Michigan
%| Fan-beam contributed by Yingying Zhang.  Basis etc. by Samuel Matej.

if nargin == 1 && streq(mask, 'test'), Gtomo_nufft_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;
[arg.nx, arg.ny] = size(mask);
arg.nb = dim_data(1);
arg.na = dim_data(2);

% option defaults
arg.chat = 0;
arg.nthread = 1;
arg.class = 'fatrix2'; % new default 2012-03-12

arg.dx = 1;
arg.ds = 1;
arg.strip_width = 1;
arg.orbit = 360;
arg.orbit_start = 0;
arg.xscale = 1;
arg.yscale = 1;

arg.dis_src_det = [];
arg.dis_iso_det = [];
arg.source_offset = 0;
arg.offset_s = 0;

% fix: we could over-sample in the future, especially if nb not a power of 2!
%arg.Krho = 2*ceil(arg.nb/2);
arg.rmax = []; % maximum object radius in FOV, relevant to sampling
arg.Krho = arg.nb;
arg.interp = {};
arg.is.kaiser = false;
arg.J = 6;
arg.Kd = 2 * size(mask);

arg.is.dsft1 = false;
arg.is.dsft2 = false;
arg.is.complex = false;
arg.is.nonneg = false; % enforce nonnegativity of output if input is nonnegative?
arg.is.shift0 = true;
arg.is.limit_pi = false;

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
arg.basis.diam	= 1;	% diameter relative to dx
arg.basis.shape	= 0;
arg.basis.m	= 0;
arg.basis.dim	= 2;
arg.basis.kernel = [];	% JxJ matrix

% model for beam radial PSF
arg.beam.type	= 'rect';	% 'rect', 'line', 'Gauss', 'KB'
arg.beam.diam	= 1;		% width of Gauss and KB relative to ds
arg.beam.shape	= 0;
arg.beam.m	= 0;

% options specified by name/value pairs
subs = {'ray_spacing', 'ds'; 'pixel_size', 'dx'; 'channel_offset', 'offset_s'};
arg = vararg_pair(arg, varargin, 'subs', subs);

% initialize geometry stuff, parallel or fan
arg = Gtomo_nufft_setup(arg);

arg.nd = prod(dim_data);
arg.np = sum(mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

% build Fatrix object
switch arg.class
case 'Fatrix'
	ob = Fatrix(dim, arg, 'caller', 'Gtomo_nufft', ...
		'forw', @Gtomo_nufft_forw, 'back', @Gtomo_nufft_back);
case 'fatrix2'
	forw = @(arg, x) Gtomo_nufft_forw_1(arg, x);
	back = @(arg, y) Gtomo_nufft_back_1(arg, y);
	warn 'todo: untested'
	ob = fatrix2('imask', mask, 'arg', arg, 'caller', 'Gtomo_nufft', ...
		'odim', dim_data, 'forw', forw, 'back', back);
otherwise
	fail('class %s unknown', arg.class)
end


% Gtomo_nufft_setup()
function arg = Gtomo_nufft_setup(arg)

arg.is.fan = ~isempty(arg.dis_src_det);
if arg.is.fan % trick: isinf([]) returns []
	arg.is.fan = arg.is.fan & ~isinf(arg.dis_src_det);
end

if arg.is.fan
	arg = Gtomo_nufft_init_fan(arg);
else
	arg = Gtomo_nufft_init_par(arg);
end

% fix: need "center_x and center_y" here, and consider flip_y or yscale
arg.nxy_shift = size(arg.mask)/2 - arg.is.shift0 * 0.5;	% shifts

%
% create the structure needed by 2D NUFFT
%
if arg.is.dsft2
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

% make filter to be used in sinogram spectral domain
if ~arg.is.fan
	arg.tomo_filter = Gtomo_nufft_filter(arg.omega, arg, ~arg.is.fan);
else
	% arg.omega is for (\rho > 0): no DC
	if mod(arg.na,2) == 0
		arg.tomo_filter = Gtomo_nufft_filter([arg.omega; -arg.omega], arg, ~arg.is.fan);
	elseif mod(arg.na,2) == 1
		arg.tomo_filter = Gtomo_nufft_filter(arg.omega, arg, ~arg.is.fan);
	end
	arg.tomo_filter = arg.tomo_filter * arg.fan.del_rho;
end

% trick: set tomo filter to pass only the "aliasing free" part
% added on 2005-8-24, because aliasing is visible in MTF
% this is a bit inefficient because we compute spectrum but zero it!
if arg.is.limit_pi
%	good = sqrt(arg.omega(:,1).^2 + arg.omega(:,2).^2) <= pi;
	good = (abs(arg.omega(:,1)) < pi) & (abs(arg.omega(:,2)) < pi);
	good = reshape(good, size(arg.omega,1)./arg.na, arg.na);
	arg.tomo_filter = arg.tomo_filter .* good;
	im(arg.tomo_filter), prompt
end


% Gtomo_nufft_init_par()
function arg = Gtomo_nufft_init_par(arg)

om_ang = ([0:(arg.na-1)]/arg.na * arg.orbit + arg.orbit_start) * pi / 180;

% issue: for parallel ray case, we limit ourselves to a FFT, which seems
% to prevent us from using the entire [-pi,pi]^2 box.
if isempty(arg.Krho)
	dx = arg.dx;
	dr = arg.ds;
	fov = dx * max(arg.nx, arg.ny);
	drho = 1 / fov; % nominal spectral sampling to avoid 2D spatial aliasing
	arg.Krho = ceil(1 / (drho * dr)); % considering 1D FFT
%	arg.Krho = ceil(sqrt(2) * arg.Krho); % trick? useless?
	arg.Krho = max(arg.Krho, arg.nb); % usually it ends up being nb
	arg.drho = 1 / (arg.Krho * dr); % hard constraint for 1D FFT
	printm('fov=%g drho=%g Krho=%d', fov, arg.drho, arg.Krho)
	om_rad = 2*pi * [-arg.Krho/2:arg.Krho/2-1] * arg.drho * dx; % per 2D NUFFT
else
	om_rad = [-arg.Krho/2:arg.Krho/2-1] * (2*pi/arg.Krho) * arg.dx / arg.ds;
end
[rr, aa] = ndgrid(om_rad, om_ang);

% desired k-space coordinates in [-pi,pi)
arg.omega = [arg.xscale * col(rr .* cos(aa)), -arg.yscale * col(rr .* sin(aa))];

% check to see if object lies within scanner FOV
fov_par_check(arg)

if arg.dx > arg.ds
	warn 'big pixels and small rays implemented inefficiently'
% fix: we could use smaller Krho then zero pad instead of tomo_filter trick
end


% Gtomo_nufft_init_fan()
function arg = Gtomo_nufft_init_fan(arg)

if isempty(arg.Krho)
	arg.Krho = arg.nb;
end

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

fan.del_sigma = arg.ds / arg.dis_src_det; % in radians
fan.sigma = ([-(arg.nb-1)/2:(arg.nb-1)/2]' ...
			- arg.offset_s) * fan.del_sigma;

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
	arg.rmax = max(abs(fan.sin_rad(:)));
end
arg.fan.del_rho = 1/(2*arg.rmax); % it has to be del_rho <= 1/dmax
% arg.fan.del_rho = 1/(2*arg.rmax*sqrt(2));
arg.fan.omegam = (-1) * 2*pi * arg.fan.del_rho * fan.sin_rad;

% center shift shouldn't be 0 since it wasn't absorbed into tomo_filter
int1c = arg.interp; % copy for 1d
if ~isempty(arg.interp) && streq(arg.interp{1}, 'table')
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
if mod(arg.Krho, 2) == 1
	k_rad = 2*pi * [-(arg.Krho - 1)/2:(arg.Krho - 1)/2]' * arg.fan.del_rho * arg.dx;
elseif mod(arg.Krho, 2) == 0
	k_rad = 2*pi * [-arg.Krho/2:arg.Krho/2]' * arg.fan.del_rho * arg.dx;
else
	warning 'mod issue'
	keyboard
end

arg.Kr = size(k_rad, 1);

if mod(arg.na, 2) == 1
	k_ang = fan.beta; % \theta when sigma=0
elseif mod(arg.na,2) == 0
	k_ang = fan.beta(1:arg.na/2); % \theta when sigma=0
	% for filter
	%[rr1 aa1] = ndgrid(k_rad((Kr - 1)/2 + 2:end), fan.beta);
	%arg.fil_omega = [arg.xscale * col(rr1 .* cos(aa1)), -arg.yscale * col(rr1 .* sin(aa1))];
end

if arg.chat
	printf('k_rad size %g %g k_ang size %g %g', ...
		size(k_rad), size(k_ang))
end

[rr, aa] = ndgrid(k_rad((arg.Kr - 1)/2 + 2:end), k_ang);

% desired k-space coordinates in [-pi,pi)
% arg.allomega = [arg.xscale * col(rr1 .* cos(aa1)), -arg.yscale * col(rr1 .* sin(aa1))];
arg.omega = [arg.xscale * col(rr .* cos(aa)), -arg.yscale * col(rr .* sin(aa))];

if arg.dx > arg.ds
	warning 'big pixels and small rays may require oversampling k'
end

% 1D NUFFT structure
if arg.is.dsft1
	% no need for 1D nufft structure
	arg.fan.n_shift = arg.Krho / 2;
else
	arg.fan.r_st = ...
	nufft_init(arg.fan.omegam, (arg.Kr-1)/2, arg.J, arg.Kr-1, -1, int1c{:}); % shift acct for index from 0, not 1
end


% fov_par_check()
function fov_par_check(arg)

x = ([0:arg.nx-1] - (arg.nx-1)/2) * arg.dx;
y = ([0:arg.ny-1] - (arg.ny-1)/2) * arg.dx;
[x, y] = ndgrid(x, y);
r = sqrt(x.^2 + y.^2);
mask_max = max(r(arg.mask(:))); % check within mask only!
r_max = ((arg.nb-1)/2 - abs(arg.offset_s)) * arg.ds; % scanner FOV radius
if mask_max > r_max
	printf('Scanner FOV radius = %g', r_max)
	printf('Object support maximum radius = %g', mask_max)
	printf 'Object support exceeds FOV of this system'
	warning 'Truncated projections probably will not work, as implemented'
end

if 0 && max(arg.nx,arg.ny) * arg.dx > arg.nb * arg.ds % old way
%	warning 'truncated projections probably will not work now'
	error 'truncated projections probably will not work now'
	% fix: could be made to work by choosing large enough Krho i think
end


% fov_fan_check()
function fov_fan_check(arg)
x = ([0:arg.nx-1] - (arg.nx-1)/2) * arg.dx;
y = ([0:arg.ny-1] - (arg.ny-1)/2) * arg.dx;
[x, y] = ndgrid(x, y);
r = sqrt(x.^2 + y.^2);
%fov_max = sqrt(nx^2 + ny^2) / 2 * arg.dx;
mask_max = max(r(arg.mask(:)));
s_max = ((arg.nb-1)/2 - abs(arg.offset_s)) * arg.ds;
r_max = arg.dis_src_iso * sin(s_max / arg.dis_src_det); % scanner FOV radius
if mask_max > r_max
	printf('Scanner FOV radius = %g', r_max)
	printf('Object support maximum radius = %g', mask_max)
	printf 'Object support exceeds FOV of this system'
	warning 'Truncated projections probably will not work, as implemented'
end


% Gtomo_nufft_init_basis()
% initialization of the image basis function (blob/Gauss/...)
% by S. Matej
function ob = Gtomo_nufft_init_basis(ob)

if streq(ob.basis.type, 'KB')
	dim1 = floor(ob.basis.diam/2);
	dim2 = 1 + 2*dim1;
	row2 = [-dim1:dim1].^2;
	onev(1:dim2) = 1;
	kernel = row2'*onev + onev'*row2;
	%normft = ...
	%kaiser_bessel_ft(0, ob.basis.diam, ob.basis.shape, ob.basis.m, ob.basis.dim)
	if ob.basis.dim==2
		kernel = sqrt(kernel);
		ob.basis.kernel = kaiser_bessel( ...
			kernel, ob.basis.diam, ob.basis.shape, ob.basis.m);

	elseif ob.basis.dim==3
		for i = 1:dim2
			kernel3(:,:,i) = kernel+row2(i);
		end
		kernel3 = sqrt(kernel3);
		ob.basis.kernel = kaiser_bessel(...
			kernel3, ob.basis.diam, ob.basis.shape, ob.basis.m);

	else
		fail('basis function dimension %g not supported', ob.basis.dim)
	end

	norm = sum(sum(sum(ob.basis.kernel)));
	ob.basis.kernel = ob.basis.kernel/norm;

elseif streq(ob.basis.type, 'Gauss')
	error('Gaussian basis not yet implemented - kernel not calculated')

elseif isempty(ob.basis.type) || streq(ob.basis.type, 'pixel') || ...
	streq(ob.basis.type, 'no')

else
	fail('basis function %s not implemented', ob.basis.type)

end


% Gtomo_nufft_forw(): y = G * x
function y = Gtomo_nufft_forw(arg, x)

[x, ei] = embed_in(x, arg.mask, arg.np);

for ll=1:size(x,3)
	y(:,:,ll) = Gtomo_nufft_forw_1(arg, x(:,:,ll));
end

y = ei.shape(y);


% Gtomo_nufft_forw_1(): y = G * x
% for a single image / sinogram
function y = Gtomo_nufft_forw_1(arg, x)

% full dsft (doesn't have nxy_shift implemented, built into tomo_filter)
if arg.is.dsft2
	y = jf_mex('dtft,forward', arg.omega', double(x), int32(arg.nthread));
	y = cast(y, class(x));

else % nufft
	y = nufft(x, arg.st);
	if arg.is.fan && mod(arg.na, 2) == 0 % fix: this logic seems wrong!
		y = [y; conj(y)]; % use half angle (0, pi) to replicate the other half (pi, 2pi)
	end
end

% for DC
if arg.is.fan
	y0 = sum(x(:)); % dc
	y0 = y0 .* Gtomo_nufft_filter([0 0], arg, ~arg.is.fan) .* arg.fan.del_rho;
end

y = Gtomo_nufft_spectral_filter(arg, y);

% parallel beam case: filter, then inverse DFT
if ~arg.is.fan % parallel

	K = arg.Krho;
%	y = y([(K/2+1):K, 1:(K/2)], :); % fft1shift
	y = fftshift(y, 1);		% fft1shift
	y = ifft(y);			% [K na] inverse 1D FFTs
	if K > arg.nb
		y = y([1:arg.nb]+[K-arg.nb]/2,:);
		warning 'over-sampling radially not fully tested'
	end


% fan beam
else
	% using 1D nufft to compute 1D nuifft
	if arg.is.dsft1
		y = jf_mex('dtft,forward', arg.fan.omegam', y, int32(arg.nthread));
		% phase effect due to n_shift:
		y = y .* repmat(exp(1i*arg.fan.omegam*arg.fan.n_shift), [1 ncol(y)]);
	else
		y = nufft(y, arg.fan.r_st); % from Krho/Krho+1 to nb
		y = 2 .* real(y) + ones(size(y)) .* y0;
	end

	% interp1 along beta to get values at desired beta locations
	y = fractional_delay(y.', -arg.delay1).'; % delay for each column

	% use linear interp to see if there is circular err pattern
%	for ib = 1:arg.nb
%		y(ib,:) = perodic_shift(y(ib,:), arg.delay1(ib));
%	end

end

if ~arg.is.complex % real
	y = real(y);
	if arg.is.nonneg && all(x(:) >= 0) % input is nonnegative
		y = max(y, 0);
	end
end


% Gtomo_nufft_back(): x = G' * y
% full backprojection
function x = Gtomo_nufft_back(arg, y)

[y, eo] = embed_out(y, [arg.nb arg.na]);

for ll=1:size(y,3)
	x(:,:,ll) = Gtomo_nufft_back_1(arg, y(:,:,ll));
end

x = eo.shape(x, arg.mask, arg.np);


% Gtomo_nufft_back_1(): x = G' * y
% for a single image / sinogram
function x = Gtomo_nufft_back_1(arg, y)

K = arg.Krho;

if ~arg.is.fan % parallel
	if K > arg.nb
		npad = (K-arg.nb)/2; % fix: probably not always correct odd/even
		y = [zeros(npad, arg.na); y; zeros(npad, arg.na)];
		warning 'over-sampling radially not fully tested'
	end
	y = (1/K) * fft(y);			% adjoint of ifft()
%	x = x([(K/2+1):K, 1:(K/2)], :);		% fft1shift
	y = fftshift(y, 1); % fft1shift

else % fan
	if K > arg.nb
%		y = [y; zeros(K-arg.nb, arg.na)];
		warning 'over-sampling radially not tested and probably wrong!'
	end

	% take fractional_delay for each sigma along beta, back to
	% uniform \theta (values at uniform \theta grid)
	y = fractional_delay(y.', arg.delay1).'; % delay for each column
	y0 = sum(y(:)) .* conj(Gtomo_nufft_filter([0 0], arg, ~arg.is.fan) .* arg.fan.del_rho);

%	for ib = 1:arg.nb
%		y(ib,:) = perodic_shift(y(ib,:), -arg.delay1(ib));
%	end

	% adjoint of 1D NUFFT along each col
	if arg.is.dsft1
		% phase effect due to n_shift:
		y = y .* repmat(exp(-1i*arg.fan.omegam*arg.fan.n_shift), [1 ncol(y)]);
		printm 'Warn: this dsft case may not be working'
		y = jf_mex('dtft,adjoint', arg.fan.omegam', y, ...
...%			int32(arg.Krho), ...
			int32((arg.Kr-1) / 2), ... % trick: new sampling!?
			int32(arg.nthread));
	else
		y = nufft_adj(y, arg.fan.r_st);
	end
end

% trick: to fix imaginary part.  subtle!
% frankly it is a bit mysterious why this works
% since "real()" is nonlinear so shouldn't have an adjoint!
if arg.is.shift0 && ~arg.is.fan
	y(1,:) = 2 * real(y(1,:));
end
y = y .* conj(arg.tomo_filter);
y = y(:);

% if mod(arg.Krho,2) == 0
%	y = y1(arg.Krho/2 + 2:end, :);
%	y0 = sum(y1(arg.Krho/2 + 1,:)); % dc
% else
%	y = y1((arg.Krho - 1)/2 + 2:end, :);
%	y0 = sum(y1((arg.Krho - 1)/2 + 1, :));
% end

% finally, back to 2D object space
if arg.is.dsft2
	if isreal(y)
		y = complex(y);
		warning 'faking sinogram complex'
	end
	x = jf_mex('dtft,adjoint', arg.omega', y, ...
			int32([arg.nx arg.ny]'), int32(arg.nthread));

else
	if arg.is.fan
		if mod(arg.na,2) == 1
			x = nufft_adj(y, arg.st);
			x = y0 + 2 * real(x);
		elseif mod(arg.na,2) == 0
			x = nufft_adj(y(1:length(y)/2)+conj(y(length(y)/2+1:end)), arg.st);
			x = y0 + 2 * real(x);
		end
	else
		x = nufft_adj(y, arg.st);
	end

	x = reshape(x, arg.nx, arg.ny);
end

if ~arg.is.complex
	x = real(x);
	if arg.is.nonneg && all(y(:) >= 0) % input is nonnegative
		x = max(x, 0);
	end
end

x = x .* arg.mask; % apply mask


% Gtomo_nufft_spectral_filter()
function y = Gtomo_nufft_spectral_filter(ob, y)

%K = ob.Krho;
%if mod(K, 2) == 1
%	y = reshape(y, [(K - 1)/2 ob.na]);
%%	y = reshape(y, K, ob.na);
%else
%keyboard
%	y = reshape(y, [], ob.na);
%%	y = reshape(y, K + 1, ob.na);
%end
try
	y = reshape(y, size(ob.tomo_filter));
catch
	printm 'size problem'
	pr size(y)
	pr size(ob.tomo_filter)
	keyboard
end
y = y .* ob.tomo_filter;

% trick: fix imaginary part.  subtle!
% see technical report for explanation of "imaginary part fix"
% x real--> y should be hermitian symmetric G(s)=conj(G(-s))
if ob.is.shift0 && ~ob.is.fan
	y(1,:,:) = 2 * real(y(1,:,:));
else
%	warning 'do i really not need this real trick?'
end
