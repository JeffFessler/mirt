  function file = xray_material_file_name(name, varargin)
%|function file = xray_material_file_name(name, [options])
%| given name of an element, such as '01' or '01-hydrogen' or 'hydrogen',
%| or of a material such as 'water', determine the full file name for the
%| mass attenuation coefficient data for that material.
%|
%| in
%|	name	element or material name
%| option
%|	'dir'	root directory for files
%| out
%|	file	rull file name
%|
%| Copyright 2006-03-30, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(name, 'test'), xray_material_file_name_test, return, end

arg.dir = [];
arg = vararg_pair(arg, varargin);
if isempty(arg.dir)
	arg.dir = [path_find_dir([filesep 'ct']) filesep 'xray-mass-atten'];
end

xm_top load
zname = xm_top('get', name);
if ~isempty(zname) % element
	file = [arg.dir filesep 'element' filesep zname];

else % compound
	file = [arg.dir filesep 'compound' filesep name];
	if exist(file, 'file'), return, end

	file = [arg.dir filesep name '.dat'];
	if ~exist(file, 'file')
		fail('mass attenuation data file "%s" not found', file)
	end
end


%
% xm_top()
%
function out = xm_top(varargin)
persistent xm_list

switch varargin{1}

case 'load'
	xm_list = {};
	xm_load

case 'add'
	xm_list{str2num(varargin{2})} = varargin{3};

case 'get'
	name = varargin{2};

	if isnumeric(name) % 1
		if 1 <= name && name <= length(xm_list)
			out = sprintf('%02d-%s', name, xm_list{name});
		else
			fail('bug for %d', name)
		end

	elseif ischar(name) && length(name) <= 2 % '01' or '1'
		out = str2num(name);
		if isempty(out), fail('bad material "%s"', name), end
		out = xm_top('get', out);

	elseif ischar(name)
		for iz=1:length(xm_list)
			out = sprintf('%02d-%s', iz, xm_list{iz});
			if streq(name, out), return, end
			if streq(name, xm_list{iz}), return, end
		end
		out = '';

	else
		out = '';
	end

otherwise
	error bug
end


function xray_material_file_name_test
xray_material_file_name('2');
xray_material_file_name(2);
xray_material_file_name('helium');
xray_material_file_name('02-helium');
xray_material_file_name('water');
xray_material_file_name('lexan');
xray_material_file_name('04-beryllium');


function xm(zn, name)
xm_top('add', zn, name)

% list of all the elements by number and name
function xm_load
xm 01 hydrogen
xm 02 helium
xm 03 lithium
xm 04 beryllium
xm 05 boron
xm 06 carbon-graphite
xm 07 nitrogen
xm 08 oxygen
xm 09 fluorine
xm 10 neon
xm 11 sodium
xm 12 magnesium
xm 13 aluminum
xm 14 silicon
xm 15 phosphorus
xm 16 sulfur
xm 17 chlorine
xm 18 argon
xm 19 potassium
xm 20 calcium
xm 21 scandium
xm 22 titanium
xm 23 vanadium
xm 24 chromium
xm 25 manganese
xm 26 iron
xm 27 cobalt
xm 28 nickel
xm 29 copper
xm 30 zinc
xm 31 gallium
xm 32 germanium
xm 33 arsenic
xm 34 selenium
xm 35 bromine
xm 36 krypton
xm 37 rubidium
xm 38 strontium
xm 39 yttrium
xm 40 zirconium
xm 41 niobium
xm 42 molybdenum
xm 43 technetium
xm 44 ruthenium
xm 45 rhodium
xm 46 palladium
xm 47 silver
xm 48 cadmium
xm 49 indium
xm 50 tin
xm 51 antimony
xm 52 tellurium
xm 53 iodine
xm 54 xenon
xm 55 cesium
xm 56 barium
xm 57 lanthanum
xm 58 cerium
xm 59 praseodymium
xm 60 neodymium
xm 61 promethium
xm 62 samarium
xm 63 europium
xm 64 gadolinium
xm 65 terbium
xm 66 dysprosium
xm 67 holmium
xm 68 erbium
xm 69 thulium
xm 70 ytterbium
xm 71 lutetium
xm 72 hafnium
xm 73 tantalum
xm 74 tungsten
xm 75 rhenium
xm 76 osmium
xm 77 iridium
xm 78 platinum
xm 79 gold
xm 80 mercury
xm 81 thallium
xm 82 lead
xm 83 bismuth
xm 84 polonium
xm 85 astatine
xm 86 radon
xm 87 francium
xm 88 radium
xm 89 actinium
xm 90 thorium
xm 91 protactinium
xm 92 uranium
