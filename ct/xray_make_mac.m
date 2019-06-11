  function mac = xray_make_mac(xrs, mas)
%|function mac = xray_make_mac(xrs, mas)
%|
%| xray_make_mac()
%| interpolate mass atten coef to source energy sampling
%| in
%|	xrs	strum	X-ray spectra; see xray_read_spectra.m
%|	mas	strum	mass attenuation coefficients; see xray_read_mac.m
%| out
%|	mac.mac	[ne L]
%|	mac.bar	[M L]	effective mass attenuation coefficents near 0
%|
if nargin < 2, help(mfilename), error(mfilename), end

mac.mac = mas.mac(xrs.en); % [ne L]
mac.bar = mas.mean(xrs.en, xrs.Ide); % [M L]

end % xray_make_mac()
