 function [mass_atten, kev, mtype, file] = ...
		xray_read_atten(mtype, kev_in, varargin)
%function [mass_atten, kev, mtype, file] = ...
%|		xray_read_atten(mtype, kev_in, [options])
%|
%| Read mass attenuation coefficients for a given material type.
%| Optionally interpolate onto desired energies.
%| 
%| in
%|	mtype			'aluminum', 'copper', 2, '2', '02-helium', ...
%|				See xray_material_file_name.m
%|		(Optionally mtype can be a cell array of several materials.
%|		If so, kev_in is mandatory and the output will be an array.)
%|	kev_in		[N 1]	optional desired energies [in keV]
%|
%| options
%|	'units'		cm | mm	default: cm
%|	'interp'	{}	interpolator type.  default {'pchip', 'extrap'}
%|	shortfile	0|1	return short file name instead of full path
%|
%| out
%|	mass_atten	[L 1]	mass attenuation coefficients [cm^2/g],
%|				If "energy" input is provided, then [N L].
%|	kev		[N 1]	at these corresponding energies.
%|	mtype		{L}	material type
%|	file		{L}	file name(s) for each material
%|
%| Copyright 2004-05-1, Jeff Fessler, University of Michigan

% default is to show example
if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(mtype, 'test'), xray_read_atten_test, return, end
if ~isvar('kev_in'), kev_in = []; end

if isnumeric(mtype) && length(mtype) > 1
	mtype = num2cell(mtype);
end

if iscell(mtype)
	if isempty(kev_in), error 'kev_in required', end
	mass_atten = zeros(length(kev_in), length(mtype)); % array
	for ll=1:length(mtype)
		[mass_atten(:,ll) dum mtype{ll} file{ll}] = ...
			xray_read_atten(mtype{ll}, kev_in, varargin{:});
	end
	kev = kev_in;
return
end

% defaults
arg.units = 'cm';
%arg.interp = {'linear', 'extrap'};
%arg.interp = {'spline', 'extrap'};
arg.interp = {'pchip', 'extrap'};
arg.shortfile = false;

arg = vararg_pair(arg, varargin);

is_mm = 0;
if streq(arg.units, 'mm')
	is_mm = 1;
elseif ~streq(arg.units, 'cm')
	error 'bad units'
end

file = xray_material_file_name(mtype);
tmp = load_ascii_skip_header(file); % read uncommented lines
kev = tmp(:,1) * 1000;	% keV
mass_atten = tmp(:,2);	% mass atten. coeff. [cm^2/g]

if arg.shortfile
	file = regexprep(file, '.*\/', '');
end

if 0 % look at k-edges
	jump = find(diff(kev) == 0);
	if ~isempty(jump)
		file = regexprep(file, '.*\/', '');
		jump = num2str(kev(jump)', '%5.1f');
		printf('%s %s', file, jump)
	end
end

% if input energies are specified, then interpolate onto those.
% trick: allow for the k-edge jumps!
if isvar('kev_in') && ~isempty(kev_in)
	mass_atten = log(mass_atten); % interpolate on a log scale
	mass_atten = interp1_jump(kev, mass_atten, kev_in, ...
		'monodown', arg.interp{:});
	mass_atten = exp(mass_atten);
	kev = kev_in;
end

if is_mm
	mass_atten = mass_atten * 100;
end


% xray_read_atten_test()
% example usage
function xray_read_atten_test
mtypes = {'lead', 'aluminum', 'water', 'lexan'};
mtypes = {'water', 'bone'};
mtypes = {'iodine', 'cesium', 'csi'};
arg1 = {};
kev = [10:510]';
for ll=1:length(mtypes)
	mtype = mtypes{ll};
	[mac1 kev1] = xray_read_atten(mtype);
	ie = min(kev) <= kev1 & kev1 <= max(kev);
	arg1 = {arg1{:}, kev1(ie), mac1(ie), 'o'};
end
mac2 = xray_read_atten(mtypes, kev); % test multiple types (cell)
if im
	clf, semilogy(arg1{:})
	legend(mtypes{:})
	hold on, semilogy(kev, mac2, '-'), hold off
	xlabel 'KeV', ylabel 'mass atten. coef.', axis tight
end
