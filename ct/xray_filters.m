  function atten = xray_filters(mtype, thickness, energy, varargin)
%|function atten = xray_filters(mtype, thickness, energy, [options])
%|
%| Compute X-ray photon survival probability as a function of energy
%| for various materials.
%| in
%|	mtype			'aluminum', 'copper', ...
%|				can be a cell array {L} for multiple filters
%|	thickness		in cm (can be an array [L] if mtype is cell)
%|	energy	[N,1]		in keV (vector)
%| option
%|	'units'	cm | mm		default: cm
%| out
%|	atten	[N,L]		unitless survival probabilities
%|
%| Copyright 2001-04-27, Jeff Fessler, University of Michigan

if nargin == 1 && streq(mtype, 'test'), xray_filters_test, return, end
if nargin < 3, ir_usage, end

arg.units = 'cm';
arg = vararg_pair(arg, varargin);

if iscell(mtype)
	LL = length(mtype);
	if length(thickness) ~= LL, error 'thickness / material mismatch', end
	atten = zeros(length(energy), LL);
	for ll=1:LL
		atten(:,ll) = xray_filters(mtype{ll}, thickness(ll), ...
			energy(:), 'units', arg.units);
	end
return
end

mass_atten = xray_read_atten(mtype, energy, 'units', arg.units); % [N,L]
density = xray_read_dens(mtype, 'units', arg.units); % [L,1]
atten = exp(-mass_atten .* thickness .* density);


function xray_filters_test
kev = [20:200]';
mtype = {'lead', 'copper', 'aluminum'};
t = xray_filters(mtype, [0.01 0.2 0.1], kev);
if im
	clf, semilogy(kev, t, '-o'), axis([20 200 10^-2 1]),
	ir_legend(mtype)
end
