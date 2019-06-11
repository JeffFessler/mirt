 function mass = xray_read_atomic_mass(mtype, varargin)
%function mass = xray_read_atomic_mass(mtype)
%|
%| Read atomic mass for a given material type.
%| 
%| in
%|	mtype			'aluminum', 'copper', 2, '2', '02-helium', ...
%|				See xray_material_file_name.m
%|		(Optionally mtype can be a cell array of several materials.
%|		If so, kev_in is mandatory and the output will be an array.)
%|
%| out
%|	mass	[M 1]	atomic mass
%|
%| Copyright 2008-04-27, Jeff Fessler, The University of Michigan

% default is to show example
if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(mtype, 'test'), xray_read_atomic_mass_test, return, end

if isnumeric(mtype) && length(mtype) > 1
	mtype = num2cell(mtype);
end

if iscell(mtype)
	mass = zeros(size(mtype)); % array
	for ll=1:length(mtype)
		mass(ll) = xray_read_atomic_mass(mtype{ll}, varargin{:});
	end
return
end

% defaults
arg.units = '';
arg = vararg_pair(arg, varargin);

% trick: mass file near is near material file, and name parsing
if isnumeric(mtype)
	file = xray_material_file_name(1); % trick to allow z > 92
	zz = mtype; % atomic number
else
	file = xray_material_file_name(mtype);

	% atomic number:
	zz = regexprep(file, '/.*\/', '');
	zz = regexprep(zz, '-.*', '');
	zz = str2num(zz);
end

% atomic mass file
file = regexprep(file, 'element.*', 'atomic-mass/atomic-mass-list');

tmp = load_ascii_skip_header(file); % read uncommented lines
if zz < 1 || zz > length(tmp)
	fail('bad element %d', zz)
end
tmp = tmp(zz,:);
if tmp(1) ~= zz, fail('bug %d', zz), end

mass = tmp(2);


% xray_read_atomic_mass_test()
% example usage
function xray_read_atomic_mass_test
mtypes = {'lead', 2, '02-helium'};
mass = xray_read_atomic_mass(mtypes);

zz = 1:118;
mass = xray_read_atomic_mass(zz);
if im
	plot(zz, zz ./ mass, '-o')
	xlabel 'Z'
	ylabel 'Z/A'
%	ir_savefig fig-z-a
end
