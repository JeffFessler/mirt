 function kspace = ir_mri_kspace_ga_radial(Nspoke, Nro, varargin)
%function kspace = ir_mri_kspace_ga_radial(Nspoke, Nro, varargin)
%|
%| Generate k-space sampling pattern for "golden angle" radial sampling.
%|
%| in
%|	Nspoke		number of spokes
%|	Nro		number of samples in each readout/spoke
%|
%| option
%|	start		first angle in series [radians] (default pi/2)
%|	angle		angular spacing [radians] (default GA)
%|	delta_ro	readout spacing (default 1/Nro)
%|	shift		shift due to gradient delays (default 0)
%|			radial sample locations are ir * delta_ro
%|			where ir = [-(Nro/2 - 1):1:Nro/2] + shift
%|	show		plot k-space locations (default false)
%|
%| out
%|	kspace	[Nro Nspoke 2] (single)
%|			kx and ky k-space locations for Nspoke*Nro samples
%|			in interval (-0.5 0.5] for default shift, delta_ro
%|			so default units are "cycles / sample"
%|
%| 2015-07 Mai Le
%| 2015-07-04 Jeff Fessler minor changes

if nargin == 1 && streq(Nspoke, 'test')
	ir_mri_kspace_ga_radial_test
	return
end
if nargin < 2, help(mfilename), error(mfilename), end

arg.shift = 0;
arg.angle = pi*(sqrt(5)-1)/2; % golden angle [radians]
arg.start = pi/2;
arg.show = false;
arg.delta_ro = 1/Nro;
arg = vararg_pair(arg, varargin);

assert(mod(Nro,2) == 0, 'Number along read out not even!')

phi = [0:Nspoke-1] * arg.angle + arg.start; % angles

rho = [-(Nro/2 - 1):1:Nro/2] * arg.delta_ro; % radial samples
rho = rho + arg.shift * arg.delta_ro; % apply gradient delay shift

phi = single(phi);
rho = single(rho);

kx = rho' * cos(phi); % [Nro Nspoke]
ky = rho' * sin(-phi);
kspace = cat(3, kx, ky); % [Nro Nspoke 2]

if arg.show && im
	plot(kx, ky, '.');
	kmax = arg.delta_ro * Nro / 2;
	axis([-1 1 -1 1] * 1.1 * kmax)
	axis square
	xtick([-1:1] * kmax), xlabelf '$k_x$'
	ytick([-1:1] * kmax), ylabelf '$k_y$'
end


% ir_mri_kspace_ga_radial_test()
% test / illustrate it, exercising most options
function ir_mri_kspace_ga_radial_test
k = ir_mri_kspace_ga_radial(13, 30, 'start', pi, ...
	'show', 1, 'shift', -0.5, 'delta_ro', 1);
