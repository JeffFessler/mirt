 function ob = spsp_Af(z, f, kz, kf, tsp, fovf, fovz, varargin)
%function ob = spsp_Af(z, f, kz, kf, tsp, fovf, fovz, varargin)
%
% input
% z : Nz by 1 vector  (cm)
% f : Nf by 1 vector  (Hz)
% kz, kf : M by 1 vector
% tsp : sampling time of RF pulse in sec
% Construct spsp_Af object
%
% This one calls the fast projections!
%
% Copyright, Sangwoo Lee, University of Michigan, 2005
% 2012-06-12, modified for fatrix2 by Jeff Fessler

arg.z = z(:);
arg.f = f(:);
arg.kz = kz(:);
arg.kf = kf(:);
arg.tsp = tsp;
arg.Nz = length(arg.z);
arg.Nf = length(arg.f);
arg.M = length(arg.kz);
arg.gamma = 26751; %rad/sec/g
arg.fovf = fovf; % Hz
arg.fovz = fovz; % cm
%arg.class = 'Fatrix'; % old
arg.class = 'fatrix2'; % new
arg.J = [8 8];
arg.oversample = 2;
arg = vararg_pair(arg, varargin);

nufft_args = {[arg.Nz arg.Nf], arg.J, arg.oversample * [arg.Nz arg.Nf], ...
	[arg.Nz arg.Nf]/2, 'table', 2^10, 'minmax:kb'};

if arg.M ~= length(arg.kf)
	error('lengths of kz and kf do not match');
end

mask = true(arg.Nz,arg.Nf);
arg.A = Gmri([kz(:)*fovz/arg.Nz, kf(:)*fovf/arg.Nf], mask, ...
	'fov', [arg.Nz arg.Nf], 'nufft', nufft_args);

switch arg.class
case 'Fatrix'
	arg.dim = [arg.Nz*arg.Nf arg.M];
	ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
		'forw', @spsp_Af_forw, 'back', @spsp_Af_back);
case 'fatrix2'
	ob = fatrix2('idim', arg.M, 'odim', arg.Nz*arg.Nf, 'arg', arg, ...
		'forw', @spsp_Af_forw, 'back', @spsp_Af_back);
otherwise
	fail('unknown class "%s"', arg.class)
end


function y = spsp_Af_forw(arg, x)
y = 1i * arg.gamma * arg.tsp * (arg.A' * x);


function x = spsp_Af_back(arg, y)
x = -1i * arg.gamma * arg.tsp * (arg.A * y);
