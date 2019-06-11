 function [f0, Bz] = mri_field_map_sim(x, z, varargin)
%function [f0, Bz] = mri_field_map_sim(x, z, varargin)
%|
%| Analytically simulate the magnetic field variation caused simple
%| geometric shapes of various susceptibilties in a uniform magnetic field.
%|
%| Currently there are two geometries impemented:
%| 1) a sphere in a uniform background based on: Yoder et al.,
%| "MRI simulator with object-specific field map calcuations,"
%| MRI, 2004.
%|
%| 2) a double-tube system based on: Lucas et al.,
%| "Diffusion imaging in the presence of static magnetic-field gradients,"
%| JMRI, 1993.
%|
%| in
%|	x	[nx]	x sample location(s) (scalar or vector)
%|	z	[nz]	z sample location(s) (scalar or vector)
%|
%| options:
%|	'type'	'sphere' (default) or 'd_tube'
%|	'gambar'	gyromagnetic ratio (default 42.576e6 Hz / T)
%|	'B0'	strength of mag. field transverse to tubes or sphere (def: 3T)
%|	'Xe'	mag. susc. of external material. (<< 1) (def: air)
%|	'Xi'	mag. susc. of internal material (sphere or ring). (<< 1) (def: water)
%|
%|	for 'sphere' (sphere)
%|		'r'	radius of sphere (def: 1)
%|		'y'	y sample location(s) (scalar or vector) (def: zeros)
%|
%|	for 'd_tube' (double tube)
%|		'ra'	inner radius of tube (def: 0.6)
%|		'rb'	outer radius of tube (def: 1.1)
%|
%| out
%|	f0	[nx ny nz]	field offset frequency over object [Hz].
%|	Bz	[nx ny nz]	field pattern over object of interest [T].
%|
%| Type mri_field_map_sim('test') for demonstration.
%|
%| The units of x, y, z, r, ra, rb must be self consistent.
%|
%| 2012-10-15 Michael Allison, University of Michigan
%| 2012-10-15 JF modified to return field map, not Bz

if nargin == 1 && streq(x, 'test'), mri_field_map_sim_test, return, end
if nargin < 2, help(mfilename), error args, end
 
% defaults
arg.type = 'sphere';
arg.B0 = 3; % Telsa
% these constants may be incorrect, says Sydney, see:
% doi 10.1016/j.mri.2003.10.001
arg.Xi = -0.721e-6; % water at 20C
arg.Xe = -9.77e-6; % air
arg.r = 1;
arg.y = zeros(size(x));
arg.ra = 0.6; % mm in paper, but unitless in eq'n
arg.rb = 1.1; % mm
arg.gambar = 42.576e6; % Hz / T

arg = vararg_pair(arg, varargin);

switch arg.type
case 'd_tube'

	r = sqrt(x.^2 + z.^2);

	Bz = (1 + arg.Xi) .* ...
		(1 - (arg.Xi./2) + ((arg.Xe - arg.Xi) .* arg.ra^2 ...
		./ (2 * r.^2)) .* (2 * z.^2 ...
		./ r.^2 - 1) ) .* arg.B0;

	% nothing specified outside of water region (ra,rb)
	Bz(find(r > arg.rb | r < arg.ra)) = arg.B0;

case 'sphere'

	Bz = (1 + (arg.Xe./3) + arg.r^3 * (arg.Xe - arg.Xi) *  ...
		(x.^2 + arg.y.^2 - 2 * z.^2) ./ ...
		(3 * (x.^2 + arg.y.^2 + z.^2).^(5/2))) .* arg.B0;

	% compute internal level with Lorentz correction.
	rad = sqrt(x.^2 + arg.y.^2 + z.^2);
	Bz(find(rad<arg.r)) = (1 + arg.Xe/3) * arg.B0;

otherwise
	fail('Unknown geometry "%s"', arg.type)
end

f0 = arg.gambar * (Bz - arg.B0); % field offset frequency in Hz

end


% mri_field_map_sim_test()
function mri_field_map_sim_test
% demonstrate the field maps using parameters from the original paper.

im plc 3 3

% test double tube
Xe = -4.8*10^-6; % Xa effective (glass and air)

x = linspace(-3,3,150);
z = linspace(-2,2,100);

[X Z] = ndgrid(x,z);

f0 = mri_field_map_sim(X, Z, 'type', 'd_tube', 'B0', 7, 'Xe', Xe);

im(1, x, z, f0, 'double tube'), cbar
axis equal, axis tight, xlabel('x'), ylabel('z')


% test sphere
B0 = 30000; %3T
x = linspace(-4,4,200);
y = linspace(-6,6,300);
z = linspace(-5,5,250);

[X Y Z] = ndgrid(x,y,z);

f0 = mri_field_map_sim(X, Z, 'type', 'sphere', 'B0', 3, 'y', Y);
%Bz = reshape(Bz,[length(y),length(x),length(z)]);

% plot central slices of sphere
ix = imin(abs(x));
iy = imin(abs(y));
iz = imin(abs(z));
im(4, x,y, f0(:,:,iz), 'sphere'), cbar, xlabel('x'), ylabel('y')
axis equal, axis tight
im(5, y,z, f0(ix,:,:), 'sphere'), cbar, xlabel('y'), ylabel('z')
axis equal, axis tight
tmp = squeeze(f0(:,iy,:))'; % trick: put z horizontal
im(6, z,x, tmp, 'sphere'), cbar, xlabel('z'), ylabel('x')
axis equal, axis tight

% plot profiles through sphere
im subplot 7, plot(squeeze(f0(:,iy,iz))), axis tight, xlabel('x')
im subplot 8, plot(squeeze(f0(ix,:,iz))), axis tight, xlabel('y')
im subplot 9, plot(squeeze(f0(ix,iy,:))), axis tight, xlabel('z')

end % mri_field_map_sim_test()
