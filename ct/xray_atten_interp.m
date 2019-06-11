 function mac = xray_atten_interp(kev, mac, kev_in, varargin)
%function mac = xray_atten_interp(kev, mac, kev_in, [options])
%|
%| Interpolate mass attenuation coefficients (mac) onto desired energies.
%| 
%| in
%|	kev	[M,1]
%|	mac	[M,1]
%|	kev_in	[N,1]		desired energies [in keV]
%|
%| option
%|	'interp' {}		default {'pchip', 'extrap'}
%| out
%|	mac	[N,1]		mass attenuation coefficients [cm^2/g],
%|
%| Copyright 2004-05-1, Jeff Fessler, University of Michigan

% default is to show example
if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(kev, 'test'), xray_atten_interp_test, return, end

% = {'linear', 'extrap'};
% = {'spline', 'extrap'};
arg.interp = {'pchip', 'extrap'};
arg = vararg_pair(arg, varargin);

% trick: allow for the k-edge jumps!
mac = log(mac); % interpolate on a log scale
mac = interp1_jump(kev, mac, kev_in, arg.interp{:});
mac = exp(mac);

% xray_atten_interp_test()
% example usage, cf Fig. 3.4 of Macovski 1983
function xray_atten_interp_test
mtype = 'water'; ax = 10.^[1 3 -2 1.];
mtype = 'lead'; ax = 10.^[1 3 -1 2.5];
[mac kev] = xray_read_atten(mtype);
kev1 = logspace(1,3,1+2^9);
mac1 = xray_atten_interp(kev, mac, kev1);
if im
	clf, loglog(kev, mac, '.', kev1, mac1, '-')
	xlabel 'KeV', ylabel 'mass attenuation coefficient [cm^2/g]'
	axis(ax)
	texts(0.7, 0.7, mtype)
	if streq(mtype, 'lead')
		text(14, 40, 'L edge', 'horizontalalignment', 'center')
		text(88, 1., 'K edge', 'horizontalalignment', 'center')
	end
end
% ir_savefig(['fig_atten_' mtype])
