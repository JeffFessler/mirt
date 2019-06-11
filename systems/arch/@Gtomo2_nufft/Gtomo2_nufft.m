 function ob = Gtomo2_nufft(nx, ny, nb, na, varargin)
%function ob = Gtomo2_nufft(nx, ny, nb, na, ...)
%
% Construct Gtomo2_nufft object, which does parallel-beam or fan-beam
% forward projection via the Fourier-slice theorem,
% with the Fourier transform computed via NUFFT.
% See Gtomo2_nufft_test.m for example usage.
% Basically, you create a system matrix object by calling:
%	G = Gtomo2_nufft( ... )
% and then you can use it thereafter by typing commands like
%	y = G * x;
% which will auto-magically perform the forward projection.
% So you can represent forward projection by a matrix multiplication
% in programs, rather than calling a subroutine.
% Why bother?  Because it is easier to debug, and more general,
% to write the routines using matrices in the first place.
%
% Besides simple utilities like display, there are the following
% capabilities of this object:
%	y = G * x		forward projection
%	x = G' * y		back projection
%
% in:
%	nb	# of radial samples in parallel beam
%		# of detector bins (sigma) in fan
%	na	# of angular samples in parallel beam
%		# of angular views (beta) in fan
%
% Optional flag arguments:
%	'fan' or 'par'	(default is parallel-beam)
%	'chat'		enable printing of debugging messages
%	'dsft'		use full DSFT rather than NUFFT
%			(this is expensive but useful for debugging)
%	'shift0'	apply shifts so center of rotation is at 0
%	'complex'	return complex projection and adjoint (unusual)
%
% Interpolation specification
%	(The default (recommended) uses minmax with KB scaling.)
%	(These others are provided for making comparisons.)
%	'kaiser'	use minmax tuned kaiser-bessel interpolator
%	'interp'	{ cell array of arguments to pass to nufft_init() }
%			(This allows complete user control of interpolator.)
%	Example: straw-man bilinear interpolation is:
%		..., 'interp', {'linear'}, 'J', 2, ...
%	Example: for straw-man uniform scaling factors in NUFFT:
%		..., 'interp', {'minmax:user', {{1}, {1}}, {{0.5}, {0.5}}}
%	Example: for table-based NUFFT for large-scale problems:
%		..., 'interp', {'table', 'minmax:kb'}
%
% Optional arguments for which a value must follow:
%	'pixel_size'	width of pixels
%	'ray_spacing'	radial sample spacing
%	'strip_width'	width of rectangular detector PSF
%	'orbit'		projection angle coverage [180 degrees]
%	'orbit_start'	first projection angle [0 degrees]
%	'xscale'	use -1 to flip in x direction
%	'yscale'	use -1 to flip in y direction
%	'mask'		binary support array
%	'J'		neighborhood size (default is 6)
%	'K1'		# of FFT points in dimension 1, default 2*nx
%	'K2'		# of FFT points in dimension 2, default 2*ny
%	'beam'		structure for beam PSF (see below)
%	'basis'		structure for image-domain basis (see below)
%	'Krho'		number of radial samples in polar space, default is nb.
%
% Required arguments for fan-beam geometry (for which a value must follow):
%	'src_det_dis'		source-to-detector distance
%	'obj2det_x'		object-to-detector distance.
%	'obj2det_y'		object-to-detector distance.
%	'source_offset'		offset between source-detector and isocenter
%	'channel_offset'	offset between central ray and detector center
%				in terms of fraction of sample(detector) spacing
%
% Copyright 2001-1, Jeff Fessler, The University of Michigan
% Extend to fan-beam geometry, 2003-11, Yingying Zhang

warning(sprintf([
'Gtomo2_nufft is obsolete and has been replaced by Gtomo_nufft\n' ...
'The arguments are *almost* the same except that now "shift0" is default.\n' ...
'And src_det_dis becomes dis_src_det and obj2det_x becomes dis_iso_det.\n' ...
'And K1 K2 becomes Kd and a mask is required and used by default.\n' ...
'And yscale no longer needs a -1 to be compatible with the usual geometry\n' ...
]))

%
% default object
%
ob = Gtomo2(0);
ob.geometry = 'par'; % default is parallel-beam

% fan-beam parameters
ob.fan.sigma = [];
ob.src_det_dis = [];
ob.obj2det_x = [];
ob.obj2det_y = [];
ob.source_offset = [];
ob.channel_offset = [];

% special attributes of this object
ob.k_ang = [];	% k-space angular samples
ob.k_rad = [];	% k-space radial samples
ob.omega = [];	% desired DTFT sample locations
ob.st = [];	% precomputed structure
ob.tomo_filter = [];	% sinogram spectrum filter
ob.Krho = [];
ob.J = 0;
ob.K1 = 0;
ob.K2 = 0;
ob.nthread = 1;			% default 1 thread for safety
ob.interp = {};			% interpolation arguments for nufft_init()
ob.is.kaiser = false;
ob.is.dsft = false;
ob.dsft = [];	% dsft matrix if requested
ob.is.complex = false;
ob.is.subref = false;	% fix: needs eliminated or replaced by index1...

% parameters related to tomographic geometry
ob.orbit = 180;			% orbit in degrees
ob.orbit_start = 0;		% first projection angle in degrees
ob.nxy_shift = [];		% image-domain shift
ob.is.shift0 = false;
ob.pixel_size = 1;
ob.ray_spacing = 1;
ob.strip_width = 1;
ob.xscale = 1;			% choose -1 to flip image left-right
ob.yscale = 1;			% choose -1 to flip image up-down

%
% image basis function used in the iterative reconstruction model (S. Matej)
% types:
%	'pixel'
%	'KB'	blob: diam=J, shape=alpha, m=KB_order, dim=dimension (2D/3D),
%		kernel - blob values (image=kernel*img_coefs, * is convolution)
%	'Gauss'	Gaussian - shape=FWHM (diam=diameter of truncation)
%	'no'	no basis function modeled/considered in the nufft
%
ob.basis.type	= 'pixel';
ob.basis.diam	= 1;	% Relative_diameter - relative to ob.pixel_size
ob.basis.shape	= 0;
ob.basis.m	= 0;
ob.basis.dim	= 2;
ob.basis.kernel	= [];	% JxJ matrix

% model for beam radial PSF
ob.beam.type	= 'rect';	% 'rect', 'line', 'Gauss', 'KB'
ob.beam.diam	= 1;		% width of Gauss and KB relative to ray_spacing
ob.beam.shape	= 0;
ob.beam.m	= 0;


if nargin < 4
	warning 'Gtomo2_nufft called with too few arguments'
	help Gtomo2_nufft
	ob = class(ob, 'Gtomo2_nufft');
	return
end

%
% required input arguments
%
ob.nx = nx;
ob.ny = ny;
ob.nb = nb;
ob.na = na;
ob.dims = [ob.nb*ob.na, ob.nx*ob.ny];

%
% default
%
ob.mask = ones(nx,ny);
ob.J = 6;
ob.K1 = 2*nx;
ob.K2 = 2*ny;

%
% optional arguments
%
if ~isempty(varargin)
	while length(varargin)
		arg = varargin{1};
		varargin = {varargin{2:end}};

		if isempty(arg)
			% do nothing for empty arguments

		elseif ~ischar(arg)
			error 'unknown non-string argument?'

		elseif streq(arg, 'fan')
			ob.geometry = 'fan';

		elseif streq(arg, 'chat')
			ob.chat = true;

		elseif streq(arg, 'complex')
			ob.is.complex = true;

		elseif streq(arg, 'kaiser')
			ob.is.kaiser = true;

		elseif streq(arg, 'dsft')
			ob.is.dsft = true;

		elseif streq(arg, 'shift0')
			ob.is.shift0 = true;

		elseif streq(arg, 'basis')
			[ob varargin] = get_arg(ob, 'basis', varargin{:});

		elseif streq(arg, 'beam')
			[ob varargin] = get_arg(ob, 'beam', varargin{:});

		elseif streq(arg, 'interp')
			[ob varargin] = get_arg(ob, 'interp', varargin{:});

		elseif streq(arg, 'mask')
			[ob varargin] = get_arg(ob, 'mask', varargin{:});
			if any(size(ob.mask) ~= [nx ny]), error 'mask size', end

		elseif streq(arg, 'orbit')
			[ob varargin] = get_arg(ob, 'orbit', varargin{:});

		elseif streq(arg, 'orbit_start')
			[ob varargin] = get_arg(ob, 'orbit_start', varargin{:});

		elseif streq(arg, 'pixel_size')
			[ob varargin] = get_arg(ob, 'pixel_size', varargin{:});

		elseif streq(arg, 'ray_spacing')
			[ob varargin] = get_arg(ob, 'ray_spacing', varargin{:});

		elseif streq(arg, 'strip_width')
			[ob varargin] = get_arg(ob, 'strip_width', varargin{:});

		elseif streq(arg, 'xscale')
			[ob varargin] = get_arg(ob, 'xscale', varargin{:});

		elseif streq(arg, 'yscale')
			[ob varargin] = get_arg(ob, 'yscale', varargin{:});

		elseif streq(arg, 'src_det_dis')
			[ob varargin] = get_arg(ob, 'src_det_dis', varargin{:});

		elseif streq(arg, 'obj2det_x')
			[ob varargin] = get_arg(ob, 'obj2det_x', varargin{:});

		elseif streq(arg, 'obj2det_y')
			[ob varargin] = get_arg(ob, 'obj2det_y', varargin{:});

		elseif streq(arg, 'source_offset')
			[ob varargin] = get_arg(ob, 'source_offset', varargin{:});

		elseif streq(arg, 'channel_offset')
			[ob varargin] = get_arg(ob, 'channel_offset', varargin{:});

		elseif streq(arg, 'J')
			[ob varargin] = get_arg(ob, 'J', varargin{:});

		elseif streq(arg, 'K1')
			[ob varargin] = get_arg(ob, 'K1', varargin{:});

		elseif streq(arg, 'K2')
			[ob varargin] = get_arg(ob, 'K2', varargin{:});

		elseif streq(arg, 'Krho')
			[ob varargin] = get_arg(ob, 'Krho', varargin{:});

		else
			error(sprintf('unknown argument <%s>', arg))
		end
	end
end

%
% polar k-space sample locations
%
if streq(ob.geometry, 'par')
	ob.k_ang = ([0:(na-1)]/na * ob.orbit + ob.orbit_start) * pi / 180;
	if isempty(ob.Krho)
		ob.Krho = nb;	% we could over-sample in the future...
				% especially if nb not a power of 2!
	end
	ob.k_rad = [-ob.Krho/2:ob.Krho/2-1] * (2*pi/ob.Krho) * ob.pixel_size / ob.ray_spacing;
	[rr aa] = ndgrid(ob.k_rad, ob.k_ang);
	% desired k-space coordinates in [-pi,pi)
	omega = [col(rr .* cos(aa)) * ob.xscale col(rr .* sin(aa)) * ob.yscale];
	% fix: should check withing mask only!
	if max(nx,ny) * ob.pixel_size > nb * ob.ray_spacing
		error 'truncated projections probably will not work now'
	end
	if ob.pixel_size > ob.ray_spacing
		warning 'big pixels and small rays may not work without oversampling k'
	end

	% fix: need "center_x and center_y" here, and consider flip_y or yscale
	ob.nxy_shift = [nx; ny]/2 - ob.is.shift0 * 0.5;	% [2,1] shifts


%
% fan beam
%
elseif streq(ob.geometry, 'fan')

	if ob.orbit ~= 360
		error 'only 360 orbit implemented due to periodic fractional delay'
	end

	if isempty(ob.Krho)
		ob.Krho = ob.nb;
	end

	%-------------------------------------------------------
	% finding the desired locations in fan(ob.na, ob.nb) data
	%-------------------------------------------------------
	ob.fan.beta = ([0:ob.na-1]'/ob.na * ob.orbit ...
			+ ob.orbit_start) * pi / 180;
	ob.fan.del_beta = diff(ob.fan.beta([1 2]));

	ob.fan.del_sigma = ob.ray_spacing / ob.src_det_dis; % in radians
	ob.fan.sigma = ([-(ob.nb-1)/2:(ob.nb-1)/2]' ...
			- ob.channel_offset) * ob.fan.del_sigma;

	%-------------------------------------------------------------------
	% find desired location in R(corresponding to sigma) for 1D NUIFFT
	% (@each angle: non in rad)
	%-------------------------------------------------------------------
	% decide from relation in PENG paper
	ob.fan.dis_src_iso = ob.src_det_dis - ob.obj2det_x;
	ob.fan.sin_rad = ob.fan.dis_src_iso * sin(ob.fan.sigma);

	r0 = abs(ob.fan.sin_rad(end) - ob.fan.sin_rad(1));
	ob.fan.del_rho = 1/r0; % it has to be del_rho <= 1/r0
	omegam = (-1) * 2*pi * ob.fan.del_rho * ob.fan.sin_rad;

	% center shift shouldn't be 0 since it wasn't absorbed into tomo_filter
	if ~isempty(ob.interp) & streq(ob.interp{1}, 'table')
		interp1 = {ob.interp{end}};
		if isnumeric(interp1{:})
			ob.fan.r_st = nufft_init(omegam, ob.Krho, ob.J, ob.Krho*2, ob.Krho/2, 'minmax:kb');
		elseif ischar(interp1{:})
			ob.fan.r_st = nufft_init(omegam, ob.Krho, ob.J, ob.Krho*2, ob.Krho/2, interp1{:});
		else
			error 'wrong interp input'
		end
	else
		ob.fan.r_st = nufft_init(omegam, ob.Krho, ob.J, ob.Krho*2, ob.Krho/2, ob.interp{:});
	end

	%------------------------------------------------------------------
	% finding corresponding k-space polar locations (uniform in both)
	%------------------------------------------------------------------
	ob.k_ang = ob.fan.beta; % \theta when sigma=0

	del_k_rad = 2*pi / ob.Krho;
	ob.k_rad = 2*pi * [-ob.Krho/2:ob.Krho/2-1]' * ob.fan.del_rho * ob.pixel_size;
	if ob.chat
		printf('k_rad size %g %g k_ang size %g %g', ...
			size(ob.k_rad), size(ob.k_ang))
	end
	[rr aa] = ndgrid(ob.k_rad, ob.k_ang);
	% desired k-space coordinates in [-pi,pi)
	omega = [	col(rr .* cos(aa)) * ob.xscale ...
			col(rr .* sin(aa)) * ob.yscale];

	% check to see if object lies within scanner FOV
	fov_fan_check(ob)

	if ob.pixel_size > ob.ray_spacing
		warning 'big pixels and small rays may require oversampling k'
	end

	% fix: need "center_x and center_y" here, and consider flip_y or yscale
	ob.nxy_shift = [nx; ny]/2 - ob.is.shift0 * 0.5;	% [2,1] shifts


else
	error 'geometry not done'
end

%
% create the structure needed by NUFFT
%

if ob.is.dsft
	ob.omega = omega;
	% no need to precompute anything
%	ob.dsft = dsft2_init(ob.omega, nx, ny, 0 !! & ob.is.shift0, ob.yscale);

else	% precompute for NUFFT

	if ob.is.kaiser
		ob.interp = {'kaiser'};
	elseif isempty(ob.interp)
		ob.interp = {'minmax:kb'};	% minmax based on KB scaling
	end

	% trick: made shift 0 since built into tomo_filter!
	if ob.chat, disp 'nufft_init', end
	ob.st = nufft_init(omega, [nx ny], [ob.J ob.J], [ob.K1 ob.K2], ...
		0 * ob.nxy_shift, ob.interp{:});

end

ob.is.empty = false;
ob = class(ob, 'Gtomo2_nufft');

if ob.chat, disp 'init_basis', end
ob = init_basis(ob);	% initialize basis function for post-reconstruction use

%
% make filter to be used in sinogram spectral domain
%
ob.tomo_filter = tomo_filter(omega, ob);

%
% optional arguments with a value, e.g., ..., 'orbit', 180, ...
%
function [ob, args] = get_arg(ob, key, varargin)
if ~length(varargin)
	error(sprintf('need value for %s', key))
end
arg = varargin{1};
ob = setfield(ob, key, arg);
args = {varargin{2:end}};
if ob.chat
	if isnumeric(arg)
		if max(size(arg)) == 1	% only print scalars
			printf('\t%s = %g', key, getfield(ob, key))
		else	% otherwise print sum (e.g., for mask)
			printf('\tsum(%s) = %g', key, sum(getfield(ob, key)))
		end
	elseif isstruct(arg)	% show first structure element if string!
		fields = fieldnames(arg);
		name = getfield(arg, fields{1});
		if ischar(name)
			printf('\t%s.%s = %s', key, fields{1}, name)
		end
	elseif ischar(arg)
		printf('\t%s = %s', key, getfield(ob, key))
	elseif isa(arg, 'cell')
		printf('\t%s {cell}', key)
	else
		warning 'display type not implemented'
	end
end


%
% fov_fan_check()
%
function fov_fan_check(ob)
x = ([0:ob.nx-1] - (ob.nx-1)/2) * ob.pixel_size;
y = ([0:ob.ny-1] - (ob.ny-1)/2) * ob.pixel_size;
[x y] = ndgrid(x, y);
r = sqrt(x.^2 + y.^2);
%fov_max = sqrt(nx^2 + ny^2) / 2 * ob.pixel_size;
mask_max = max(r(ob.mask(:)));
s_max = ((ob.nb-1)/2 - abs(ob.channel_offset)) * ob.ray_spacing;
r_max = ob.fan.dis_src_iso * sin(s_max / ob.src_det_dis); % scanner FOV radius
if mask_max > r_max
	printf('Scanner FOV radius = %g', r_max)
	printf('Object support maximum radius = %g', mask_max)
	printf 'Object support exceeds FOV of this system'
	warning 'Truncated projections probably will not work, as implemented'
end


%
% DSFT: 2d discrete-space Fourier transform matrix
% The exact form of DSFT - not sparse, so slow and big
%
function W = dsft2_init_NOT_USED(omega, nx, ny, is_shift0, yscale)

ix = ([-nx/2:nx/2-1] + is_shift0 * 0.5) * xscale;
iy = ([-ny/2:ny/2-1] + is_shift0 * 0.5) * yscale;
[ix iy] = ndgrid(ix, iy);

W = exp(-1j * (omega(:,1) * ix(:)' + omega(:,2) * iy(:)'));
