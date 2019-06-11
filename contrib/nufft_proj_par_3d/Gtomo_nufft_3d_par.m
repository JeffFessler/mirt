 function ob = Gtomo_nufft_3d_par(mask, dim_data, varargin)
%function ob = Gtomo_nufft_3d_par(mask, dim_data, options)
%
% Construct Gtomo_nufft_3d_par object. Currently only handle Fourier-based 
% forward projection for 3D parallel beam case.
% 
% See Gtomo_nufft_3d_par_test.m for example usage.
%
% Basically, you create an object calling:
%		G = Gtomo_nufft_3d_par(...)
% and then you can use it thereafter by typing commands like
%		y = G * x;
% which will auto-magically do the multiplication by calling a mex file.
%
% in 3D
%	mask	 [nx,ny,nz]	logical array of object support.
%	dim_data [4]	    data dimensions, usually: [ns nt na npo], where
%		     ns	        number of horizontal samples
%            nt         number of vertical samples
%		     na	        number of azimuthal angles 
%            npo        number of polar angles
%                       Note: now set to 1 (theta=0) for comparison with 
%                             analytical result from ellipsoid_proj.m
%
% optional arguments for which a value must follow:
%	'chat'			verbose printing of debug messages
%	'orbit_a'		azimuthal angle coverage [180 degrees]
%	'orbit_a_start'	first azimuthal angle [0 degrees]
%   'orbit_p'       polar angle coverage [-orbit_p orbit_p]
%	'xscale'	    use -1 to flip in x direction
%	'yscale'	    use -1 to flip in y direction
%   'zscale'        use -1 to flip in z direction
%   'dx'            voxel(cubic) width
%   'ds'            horizontal sample spacing 
%   'dt'            vertical sample spacing
%   'rect_s'        horizontal width of cuboid detector PSF
%   'rect_t'        vertical width of cuboid detector PSF
%	
%   'beam'		    structure for beam PSF (see below)
%	'basis'		    structure for image-domain basis (see below)
%
%	'is.dsft1'	    use full 1D DSFT rather than 1D NUFFT
%	'is.dsft2'	    use full 2D DSFT rather than 2D NUFFT
%   'is.dsft3'	    use full 2D DSFT rather than 3D NUFFT
%			(these are expensive but useful for debugging)
%	'is.shift0'	    apply shifts so center of rotation is at 0 (is default)
%	'is.complex'	return complex projection and adjoint (unusual)
%
% required arguments for 3D parallel beam geometry (for which a value must follow):
%	'offset_s' | 'channel_offset'	offset between central ray and detector center
%				in terms of fraction of sample(detector) spacing
%   'offset_t' 
%
% interpolation specification
%
%	'J'		  neighborhood size (default is 5)
%	'Kd'      # of FFT points along each dimension 2*[nx ny nz]
%	'Kv1'     number of horizontal samples in Fourier space, default is ns
%   'Kv2'     number of vertical samples in Fourier space, default is nt
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
%		..., 'interp', {'table', 2^11, 'minmax:kb'}
%
% out
%	ob [nd,np]	np = sum(mask(:)), so it is already "masked"
%
% based on Gtomo_nufft.m by Jeff Fessler and Yingying Zhang
% modified for 3D parallel beam case by Yong Long, 2007-07-27


if nargin == 1 && streq(mask, 'test'), Gtomo_nufft_3d_par_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.mask = mask;
[arg.nx arg.ny arg.nz] = size(mask);
arg.ns = dim_data(1);
arg.nt = dim_data(2);
arg.na = dim_data(3);
arg.npo = dim_data(4);

% option defaults
arg.chat = 0;
arg.nthread = 1;

arg.rect_s = 1;
arg.rect_t = 1;
arg.dx = 1;
arg.dz = 1;
arg.ds = 1;
arg.dt = 1;
arg.orbit_a = 360;
arg.orbit_a_start = 0;
arg.orbit_p = 30;
arg.xscale = 1;
arg.yscale = 1;
arg.zscale = 1;

arg.dis_src_det = [];
arg.dis_iso_det = [];
arg.source_offset = 0;
arg.offset_s = 0;
arg.offset_t = 0;

arg.Kv1 = arg.ns;
arg.Kv2 = arg.nt;

arg.interp = {};
arg.is.kaiser = logical(0);
arg.J = 5;
arg.Kd = 2 * size(mask);

arg.is.dsft1 = logical(0);
arg.is.dsft2 = logical(0);
arg.is.dsft3 = logical(0);
arg.is.complex = logical(0);
arg.is.shift0 = logical(1);
arg.is.limit_pi = logical(0);

arg.basis.type = 'cube';
arg.basis.diam	= 1;	% diameter relative to dx
arg.basis.shape	= 0;
arg.basis.m	= 0;
arg.basis.dim	= 3;
arg.basis.kernel = [];	% JxJ matrix

% model for beam radial PSF
arg.beam.type	= 'rect2';	
arg.beam.diam	= 1;		
arg.beam.shape	= 0;
arg.beam.m	= 0;

% options specified by name/value pairs
subs = {'ray_spacing', 'ds'; 'pixel_size', 'dx'; 'channel_offset', 'offset_s'};
arg = vararg_pair(arg, varargin, 'subs', subs);

% initialize geometry stuff, 3D parallel beam
arg = Gtomo_nufft_3d_par_setup(arg);

arg.nd = prod(dim_data);
arg.np = sum(mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtomo_nufft', ...
	'forw', @Gtomo_nufft_forw, 'back', @Gtomo_nufft_back);


%
% Gtomo_nufft_3d_par_setup()
%
function arg = Gtomo_nufft_3d_par_setup(arg)

arg = Gtomo_nufft_init_3d_par(arg);

arg.nxyz_shift = size(arg.mask)/2 - arg.is.shift0 * 0.5;	% shifts

%
% create the structure needed by 1D NUFFT + 2D NUFFT
%
% precompute for NUFFT
if arg.is.kaiser
    arg.interp = {'kaiser'};
elseif isempty(arg.interp)
    arg.interp = {'minmax:kb'};	% minmax based on KB scaling
end

% trick: made shift 0 since built into tomo_filter!
if arg.chat, disp 'nufft_init', end

% 1D NUFFT structure along Wz direction
arg.st.z = nufft_init(arg.omega.z, arg.nz, arg.J, arg.Kd(3), ...
    0,'minmax:kb');% arg.interp{:});

% 2D NUFFT structure on Wx-Wy plane
% only table-based NUFFT used here
% all Wx-Wy planes use the same structure, saving memory
nsa = arg.Kv1 * arg.na;
arg.st.xy = nufft_init(zeros(nsa,2), [arg.nx arg.ny], [arg.J arg.J],...
    arg.Kd(1:2), 0 * arg.nxyz_shift(1:2), arg.interp{:});

if arg.chat, disp 'init_basis', end
arg = Gtomo_nufft_init_basis(arg);	% space-domain basis function

%
% make filter to be used in sinogram spectral domain
%
arg.tomo_filter = Gtomo_nufft_filter_3D_par(arg.omega.xyz, arg, 1);

% trick: set tomo filter to pass only the "aliasing free" part
% added on 2005-8-24, because aliasing is visible in MTF
% this is a bit inefficient because we compute spectrum but zero it!
if arg.is.limit_pi
    good = (abs(arg.omega.xyz(:,1)) < pi) & (abs(arg.omega.xyz(:,2)) < pi) ...
        & (abs(arg.omega.xyz(:,3)) < pi);
    good = reshape(good, arg.Kv1, arg.na, (arg.Kv2 + 1)/2, arg.npo);
    good = permute(good, [1 3 2 4]);
    
    arg.tomo_filter = arg.tomo_filter .* good;
	im(arg.tomo_filter), prompt
end

%
% Gtomo_nufft_init_3d_par()
% 3D parllel-beam
%
function arg = Gtomo_nufft_init_3d_par(arg)
    
if isempty(arg.Kv1)
	arg.Kv1 = arg.ns;
end

if isempty(arg.Kv2)
	arg.Kv2 = arg.nt;
end

par.phi = ([0:(arg.na-1)]/arg.na * arg.orbit_a + arg.orbit_a_start) * pi / 180;
par.del_phi = diff(par.phi([1 2]));

%par.theta = ([0:(arg.npo-1)]/arg.npo * 2 * arg.orbit_p - arg.orbit_p) * pi / 180;
%par.del_theta = diff(par.theta([1 2]));
par.theta = 0; %currently only theta=0 for comparison with analytical result

arg.par.del_v1 = 1/(arg.Kv1 * arg.ds);
arg.par.del_v2 = 1/(arg.Kv2 * arg.dt);

%--------------------------------------------------------------------------
% find corresponding k-space frequency locations
%--------------------------------------------------------------------------
if mod(arg.Kv1, 2) == 1
	k_v1 = 2*pi * [-(arg.Kv1 - 1)/2:(arg.Kv1 - 1)/2]' * arg.par.del_v1 * arg.dx;
elseif mod(arg.Kv1, 2) == 0
    k_v1 = 2*pi * [-arg.Kv1/2:arg.Kv1/2]' * arg.par.del_v1 * arg.dx;
    %the number of k_v1 is arg.Kv1+1 in order to utilize hermitian symmetry
else
	warning 'mod issue'
	keyboard
end

if mod(arg.Kv2, 2) == 1
	k_v2 = 2*pi * [-(arg.Kv2 - 1)/2:(arg.Kv2 - 1)/2]' * arg.par.del_v2 * arg.dx;
elseif mod(arg.Kv2, 2) == 0
    k_v2 = 2*pi * [-arg.Kv2/2:arg.Kv2/2]' * arg.par.del_v2 * arg.dx;
    %the number of k_v2 is arg.Kv2+1 in order to utilize hermitian symmetry
else
	warning 'mod issue'
	keyboard
end

arg.Kv1 = size(k_v1, 1);
arg.Kv2 = size(k_v2, 1);

if arg.chat
	printf('k_v1 size %g %g k_v2 size %g %g', size(k_v1), size(k_v2))
end

% desired K-space coordinates in [-pi,pi)
[vv1 aa vv2 theta] = ndgrid(k_v1, par.phi, ...
    k_v2((arg.Kv2 - 1)/2 + 1:end), par.theta);
clear k_v1
arg.omega.xyz = zeros(arg.Kv1 * (arg.Kv2+1)/2 * arg.na * arg.npo, 3);
arg.omega.xyz(:,1) = arg.xscale * col(vv1 .* cos(aa) + vv2 .* sin(aa) .* sin(theta));
arg.omega.xyz(:,2) = -arg.yscale * col(vv1 .* sin(aa) - vv2 .* cos(aa) .* sin(theta));
clear vv1 aa
arg.omega.xyz(:,3) = arg.zscale * col(vv2 .* cos(theta));
clear vv2 theta

%desired Wz coordinates in K-space
[vv2 theta] =ndgrid(k_v2((arg.Kv2 - 1)/2 + 1:end), par.theta);
arg.omega.z = arg.zscale * col(vv2 .* cos(theta));
clear vv2 theta k_v2
%--------------------------------------------------------------------------

%
% Gtomo_nufft_init_basis()
% only cube now (rect * rect * rect)

function ob = Gtomo_nufft_init_basis(ob)

if streq(ob.basis.type, 'cube')
else
	error(sprintf('basis function %s not implemented', ob.basis.type))
end


%
% Gtomo_nufft_forw(): y = G * x
%
function yproj = Gtomo_nufft_forw(arg, x)

% if needed, convert concise column to 3d array
flag_column = 0;
if size(x,1) == arg.np
	flag_column = 1;
	x = embed(x, arg.mask);
end

% nufft(1D + 2D)
cpu tic
y = nufft_3d_par(x, arg.st, arg.omega.xyz(:,1:2)); 
cpu toc 'NUFFT time:'

% arrange data to format: (v1, v2, phi, theta)
y = reshape(y, [arg.Kv1, arg.na, (arg.Kv2 + 1)/2, arg.npo]);
y = permute(y, [1 3 2 4]);

y = Gtomo_nufft_spectral_filter(arg, y); 
          
% utilizing Hermitian symmetry property to get the values at another half
% frequency locations
yproj = zeros(arg.Kv1, arg.Kv2, arg.na, arg.npo);
for jj = 1 : arg.npo          
    for ii = 1 : arg.na
        y1 = flipud(fliplr(conj(y(:,2:end,ii,jj))));
        yproj(:,:,ii,jj) = [y1 y(:,:,ii,jj)];
    end
end
y = yproj(1:end-1,1:end-1,:,:); %eliminate redundent values 
 
% inverse 2D FFT, for each combination of phi and theta
yproj = zeros(arg.ns, arg.nt, arg.na, arg.npo);
for jj = 1 : arg.npo
    for ii = 1 : arg.na
        y1 = y(:,:,ii,jj);
        y1 = fftshift(y1);
        y1 = ifft(ifft(y1, [], 1), [], 2);
        yproj(:,:,ii,jj) = y1;
    end
end
 

if ~arg.is.complex
	yproj = real(yproj);
end

iz = yproj < 0;
yproj(iz) = 0;

if flag_column % column in yields column out.
    yproj = reshape(yproj, [arg.ns arg.nt arg.na arg.npo]);
end

%
% Gtomo_nufft_spectral_filter()
%
function y = Gtomo_nufft_spectral_filter(ob, y)

try
y = reshape(y, size(ob.tomo_filter));
catch
	size(y), size(ob.tomo_filter)
	keyboard
end
y = y .* ob.tomo_filter;
