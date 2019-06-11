 function [kspace, omega, wi] = mri_trajectory(ktype, arg_traj, N, fov, arg_wi)
%function [kspace, omega, wi] = mri_trajectory(ktype, arg_traj, N, fov, arg_wi)
%|
%| generate kspace trajectory samples and density compensation functions.
%|
%| in
%|	ktype		string	k-space trajectory type.  see choices below.
%|	arg_traj	cell	arguments for a specific trajectory
%|	N	[1 2|3]		target image size
%|	fov	[1 2|3]		field of view in x and y (and z)
%|	arg_wi		cell	options to pass to ir_mri_density_comp
%|
%| out
%|	kspace	[Nk 2|3]	kspace samples in units 1/fov
%|	omega	[Nk 2|3]	trajectory samples over [-pi,pi]
%|	wi	[Nk 1]		(optional) density compensation factors
%|
%| trajectory types:
%| 'cartesian' 'radial' 'cart:y/2' 'random'
%| 'half+8' 'epi-sin'
%| 'spiral0' 'spiral1' 'spiral3'
%| 'rosette3'
%| 'epi-under'
%| 'gads' % emulate golden-angle data sharing per winkelmann:07:aor
%|
%| Copyright 2004-4-21, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ktype, 'test'), mri_trajectory_test, return, end
if nargin < 4, help(mfilename), fail 'args', end
if ~isvar('arg_wi') arg_wi = {}; end

if length(N) == 1, N = [N N]; end
if length(fov) == 1, fov = fov * ones(size(N)); end


% default density compensation, which works for ordinary Cartesian (only)
wi = [];
if isempty('arg_wi')
	wi = 1 / prod(fov);
end

%
% trajectory choices
%

% ideal cartesian
switch ktype
case 'cartesian'
	if any(rem(N, 2)), warn 'odd', end % todo: antonis will fix
	o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
	o2 = ([0:(N(2)-1)]/N(2) - 0.5)*2*pi;
	if length(N) == 2
		[o1 o2] = ndgrid(o1, o2);
		omega = [o1(:) o2(:)];
	elseif length(N) == 3
		o3 = ([0:(N(3)-1)]/N(3) - 0.5)*2*pi;
		[o1, o2, o3] = ndgrid(o1, o2, o3);
		omega = [o1(:) o2(:) o3(:)];
	else
		fail 'only 2d and 3d done'
	end
	wi = 1 / prod(fov);

case 'epi-under'
	[omega, wi] = mri_trajectory_epi_under(N(1:2), fov(1:2), arg_traj{:});

case 'gads'
	if length(N) ~= 2, fail 'only 2d done', end
	[omega, wi] = mri_trajectory_gads(N, fov, arg_traj{:});

case 'radial'
	[omega, wi] = mri_trajectory_radial(N(1:2), fov(1:2), arg_traj{:});
	if length(N) == 3 % 'stack of radials'
		omega = mri_trajectory_stack(omega, N(3));
		wi = repmat(wi, N(3), 1);
		wi = wi(:) / fov(3);
	end
	if ~isempty(arg_wi)
		wi = []; % discard the default "analytical" DCFs
	end

case 'rosette3'
	if length(N) ~= 3, fail 'only 3d done', end
	[omega, wi] = mri_trajectory_rosette3(N, fov, varargin)

% half cartesian + 8 rows
case 'half+8'
	o1 = ([0:(N(1)-1)]/N(1) - 0.5)*2*pi;
	o2 = [-N(2)/2:8]/N(2) * 2*pi;
	[oo1 oo2] = ndgrid(o1, o2);
	omega = [oo1(:), oo2(:)];

% echo-planar with sinusoid:
case 'epi-sin'
	if isempty(arg_traj)
		oversample = 1;
	elseif iscell(arg_traj) && length(arg_traj) == 1
		oversample = arg_traj{1};
	else
		fail 'bad trajectory argument'
	end
	Npt = oversample*prod(N);
	t = [0:(Npt-1)]'/Npt;
	omega = [pi*sin(2*pi*t*N(2)/2) t*2*pi-pi];

% bad spiral:
case 'spiral0'
	Nspiral = 400;
	omega = linspace(0, 10*2*pi, Nspiral)';
	omega = pi*[cos(omega) sin(omega)] .* omega(:,[1 1])/max(omega);
	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
	end

% crude spiral:
case 'spiral1'
	Nspiral = round(prod(N) * pi/4);
	omega = linspace(0, N(1)*2*pi, Nspiral)';
	omega = pi*[cos(omega) sin(omega)] .* omega(:,[1 1])/max(omega);
	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
	end

% 3T spiral:
case 'spiral3'
	if fov(1) ~= fov(2) || N(1) ~= N(2), fail 'only square done', end
	[kspace, omega] = mri_kspace_spiral('N', max(N(1:2)), ...
				'fov', max(fov(1:2)));

	if length(N) == 3, % stack of these spirals
		omega = mri_trajectory_stack(omega, N(3));
	end

	if isempty(arg_wi)
		wi = abs(omega(:,1) + 1i * omega(:,2)); % simple |r| weighting
		wi = wi / (fov(1)*fov(2)); % approximate scaling
		if length(N) == 3
			wi = wi / fov(3); % cartesian-style weighting in z
		end
	end

% random
case 'random'
	rng(0)
	omega = (rand(N(1)*N(2)*2, 2)-0.5)*2*pi;

% 2D FT, undersampled in "y" (phase encode) direction
case 'cart:y/2'
	o1 = ([0:(N(1)/1-1)]/(N(1)/1) - 0.5)*2*pi;
	o2 = ([0:(N(2)/2-1)]/(N(2)/2) - 0.5)*2*pi;
	[oo1 oo2] = ndgrid(o1, o2);
	omega = [oo1(:), oo2(:)];

otherwise
	fail('unknown trajectory "%s"', ktype)
end

% convert to physical units
kspace = zeros(size(omega), 'single');
for id=1:length(N)
	dx = fov(id) / N(id);
	kspace(:,id) = omega(:,id) / (2*pi) / dx;
end

if ~isempty(arg_wi) && nargout > 2 && isempty(wi)
	wi = ir_mri_density_comp(kspace, arg_wi{:});
end


% mri_trajectory_stack()
% make 3D "stack of ..." trajectory from a 2D trajectory
function omega = mri_trajectory_stack(omega2, N3)
o3 = single([0:(N3-1)]/N3 - 0.5)*2*pi;
o3 = repmat(o3, nrow(omega2), 1); % [N12,N3]
omega = repmat(omega2, N3, 1); % [N12*N3,2]
omega = [omega o3(:)]; % [N12*N3,3]


% mri_trajectory_epi_under()
% EPI, with optional under-sampling
function [omega wi] = mri_trajectory_epi_under(N, fov, varargin)
arg.samp = true(N(2),1); % default keeps all phase-encode samples
arg = vararg_pair(arg, varargin);
nx = N(1);
ny = N(2);
omx = single([-nx/2:nx/2-1]') / nx * 2*pi; % [-pi,pi) in x
omega = [];
for iy=1:ny % # of possible phase encodes
	if arg.samp(iy)
		omy = (iy-1-ny/2) / ny * 2 * pi;
		omega = [omega; [omx omy*ones(nx,1)]];
		omx = flipud(omx); % trick: EPI sweeps back and forth
	end
end
wi = ones(size(omega,1), 1, 'single') / prod(fov);


% mri_trajectory_gads()
% emulate 2D golden angle radial sampling with data sharing
function [omega wi] = mri_trajectory_gads(N, fov, varargin)
rng(0)
arg.Nro = max(N); % # of samples in each readout/spoke
arg.delta_ro = 1 / arg.Nro;
%arg.shift = 0;
arg.shift = -0.75; % shift along read-out due to gradient delays (stress)
arg.kmax_frac = [0.20 0.35 0.501]; % fractions of maximum krad (0.5) for rings
arg.nspoke = floor(pi * arg.kmax_frac * arg.Nro);
arg.under = [1 1 0.6]; % under-sampling factor for each annulus
arg.nspoke = arg.nspoke .* arg.under;
if 1 % make fibonacci for more uniform coverage
	pr arg.nspoke
	phi = (1 + sqrt(5)) / 2;
	n = round(log(arg.nspoke * sqrt(5) + 0.5) / log(phi));
	arg.nspoke = floor(phi.^n / sqrt(5) + 0.5);
	pr arg.nspoke
end
arg.start = [0 1 2] * pi/4;
arg = vararg_pair(arg, varargin);
nring = numel(arg.nspoke);
omega = [];
kmax_frac = [0 arg.kmax_frac];
for ir=1:nring
	kspace = ir_mri_kspace_ga_radial(arg.nspoke(ir), arg.Nro, ...
		'delta_ro', arg.delta_ro, 'shift', arg.shift, ...
		'start', arg.start(ir));
	kspace = reshape(kspace, [], 2);
	krad = sqrt(sum(kspace.^2, 2));
	good = (kmax_frac(ir) <= krad) & (krad < kmax_frac(ir+1));
	kspace = kspace(good,:);
	omega = [omega; 2*pi*kspace];
end
wi = []; % no default DCF


% mri_trajectory_radial()
% todo: generalize to 3D using barger:02:trc
function [omega wi] = mri_trajectory_radial(N, fov, varargin)
arg.na_nr = 2*pi;	% default ensures proper sampling at edge of k-space
arg.na = [];		% angular spokes (default: na_nr * nr)
arg.nr = max(N)/2;	% radial samples per spoke
arg.ir = [];		% default: 0:nr
arg.omax = pi;		% maximum omega
arg = vararg_pair(arg, varargin);
if isempty(arg.ir), arg.ir = [0:arg.nr]; end
if isempty(arg.na), arg.na = 4*ceil(arg.na_nr * arg.nr/4); end % mult of 4
om = arg.ir/arg.nr * pi;
ang = [0:arg.na-1]/arg.na * 2*pi;
[om ang] = ndgrid(om, ang); % [nr+1, na]
omega = [col(om.*cos(ang)) col(om.*sin(ang))];

% density compensation factors based on "analytical" voronoi
if any(fov ~= fov(1)), fail('only square FOV implemented for radial'), end
du = 1/fov(1); % assume this radial sample spacing
wi = pi * du^2 / arg.na * 2 * arg.ir(:); % see lauzon:96:eop, joseph:98:sei
wi(arg.ir == 0) = pi * (du/2)^2 / arg.na; % area of center disk
wi = repmat(wi, [1 arg.na]);
wi = wi(:);


% mri_trajectory_rosette3()
% 3d rosette, with default parameters from bucholz:08:miw
function [omega wi] = mri_trajectory_rosette3(N, fov, varargin)
arg.f1 = 211;
arg.f2 = 117.13;
arg.f3 = 73.65;
arg.nshot = 32; % todo: shots vs arms
arg.omax = pi; % maximum omega
arg.nt = 16384; % # time samples (65.536 ms for 4 usec dt)
arg.dt = 4e-6; % # time sample spacing (4 usec)
arg.ti = []; % # time samples
arg = vararg_pair(arg, varargin);
if isempty(arg.ti)
	arg.ti = [0:arg.nt-1]' * arg.dt;
end
tmp = 2 * pi * arg.ti;
p1 = f1 * tmp;
p2 = f2 * tmp;
p3 = f3 * tmp;
kx = arg.omax * sin(p1) .* cos(p2) .* cos(p3);
ky = arg.omax * sin(p1) .* sin(p2) .* cos(p3); 
kz = arg.omax * sin(a1) .* sin(a3);
omega = [kx ky kz];
for is=1:(arg.nshot-1) % n-shot, rotate kx,ky by 2 pi / N
	ang = is * 2 * pi / arg.nshot;
	c = cos(ang);
	s = sin(ang);
	ox = c * kx + s * ky;
	oy = -s * kx + c * ky;
	omega = [omega; [ox oy kz]];
end

wi = omax^3 * abs( sin(p1)^2 .* cos(p1) .* cos(p3) ); % from bucholz:08:miw


% mri_trajectory_test
% test routine
function mri_trajectory_test

ig = image_geom_mri('nx', 2^6, 'ny', 2^6-0, 'fov', 250); % 250 mm FOV
N = ig.dim;

arg_tr = {};
arg_wi = {};
%ktype = 'cartesian';
%ktype = 'spiral';
ptype = '.';
%ktype = 'epi-sin'; arg_tr = {2};
%ktype = 'radial'; arg_tr = {'na_nr', pi/2}; ptype = '.-';
ktype = 'gads'; arg_wi = {'voronoi'};
%arg_wi = {'voronoi'};
[kspace, omega, wi] = mri_trajectory(ktype, arg_tr, N, ig.fov, arg_wi);

im pl 2 2
if im
	im subplot 1, plot(omega(:,1), omega(:,2), ptype)
	titlef('"%s" with %d k-space samples', ktype, size(omega,1))
	axis_pipi, axis square
end

printm 'setup Gnufft object'
A = Gnufft(ig.mask, ...
	{omega, N, [6 6], 2*N, [N/2], 'table', 2^10, 'minmax:kb'});

printm 'setup data'
obj = mri_objects('fov', ig.fov, 'rect2half');
xt = obj.image(ig.xg, ig.yg);
xt(end/2, end/2) = 0;
yi = obj.kspace(kspace(:,1), kspace(:,2));

printm 'conj. phase reconstruction'
xcp = A' * (wi .* yi); % apply DCF for CP
xcp = ig.embed(xcp);

im(2, ig.x, ig.y, xt, 'f true'), cbar
im(3, ig.x, ig.y, abs(xcp), 'Conj. Phase Recon.'), cbar

im subplot 4
ix = 1:ig.nx; iy = ig.ny/2+1;
if im
	plot(	ig.x, xt(ix,iy), '-', ...
		ig.x, real(xcp(ix,iy)), 'g.-', ...
		ig.x, imag(xcp(ix,iy)), 'y.-')
	axis tight
end
