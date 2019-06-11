 function rho = element_density(elem)
%function rho = element_density(elem)
% density of an element (in g/cc)
if ~nargin, help(mfilename), error(mfilename), return, end

% http://environmentalchemistry.com/yogi/periodic/density.html 
% awk '{print $5, $4, $3, $1, $2}' < tt | sort -n > element_density.m

all = load;

if isnumeric(elem)
	elem = sprintf('%02d', elem); % 1 to '01'
end
if length(elem) > 3 && elem(3) == '-'
	elem = elem(1:2); % '01-hydrogen' to '01'
end

for ii=1:length(all)
	[z abbrev name rho unit] = deal(all{ii}{:});
	rho = str2num(rho);
	if streq(elem, z) || strcmpi(elem,abbrev) || strcmpi(elem,name)
		if ~streq(unit, 'g/cc')
			error('bad units "%s" for element "%s"', unit, elem)
		end
		return
	end
end
error('unknown density for %s', elem)


%
% add()
%
function out = add(varargin)

persistent all
if streq(varargin{1}, 'reset')
	all = {};
end

if nargout
	out = all;
else
	all{end+1} = varargin;
end


%
% load()
%
function all = load

add('reset')

add 01 H Hydrogen 0.0899 g/L
add 02 He Helium 0.1785 g/L
add 03 Li Lithium 0.534 g/cc
add 04 Be Beryllium 1.848 g/cc
add 05 B Boron 2.34 g/cc
add 06 C Carbon 2.26 g/cc
add 07 N Nitrogen 1.2506 g/L
add 08 O Oxygen 1.429 g/L
add 09 F Fluorine 1.696 g/L
add 10 Ne Neon 0.9 g/L
add 11 Na Sodium 0.971 g/cc
add 12 Mg Magnesium 1.738 g/cc
add 13 Al Aluminum 2.702 g/cc
add 14 Si Silicon 2.33 g/cc
add 15 P Phosphorus 1.82 g/cc
add 16 S Sulfur 2.07 g/cc
add 17 Cl Chlorine 3.214 g/L
add 18 Ar Argon 1.7824 g/L
add 19 K Potassium 0.862 g/cc
add 20 Ca Calcium 1.55 g/cc
add 21 Sc Scandium 2.99 g/cc
add 22 Ti Titanium 4.54 g/cc
add 23 V Vanadium 6.11 g/cc
add 24 Cr Chromium 7.19 g/cc
add 25 Mn Manganese 7.43 g/cc
add 26 Fe Iron 7.874 g/cc
add 27 Co Cobalt 8.9 g/cc
add 28 Ni Nickel 8.9 g/cc
add 29 Cu Copper 8.96 g/cc
add 30 Zn Zinc 7.13 g/cc
add 31 Ga Gallium 5.907 g/cc
add 32 Ge Germanium 5.323 g/cc
add 33 As Arsenic 5.72 g/cc
add 34 Se Selenium 4.79 g/cc
add 35 Br Bromine 3.119 g/cc
add 36 Kr Krypton 3.75 g/L
add 37 Rb Rubidium 1.63 g/cc
add 38 Sr Strontium 2.54 g/cc
add 39 Y Yttrium 4.47 g/cc
add 40 Zr Zirconium 6.51 g/cc
add 41 Nb Niobium 8.57 g/cc
add 42 Mo Molybdenum 10.22 g/cc
add 43 Tc Technetium 11.5 g/cc
add 44 Ru Ruthenium 12.37 g/cc
add 45 Rh Rhodium 12.41 g/cc
add 46 Pd Palladium 12.02 g/cc
add 47 Ag Silver 10.5 g/cc
add 48 Cd Cadmium 8.65 g/cc
add 49 In Indium 7.31 g/cc
add 50 Sn Tin 7.31 g/cc
add 51 Sb Antimony 6.684 g/cc
add 52 Te Tellurium 6.24 g/cc
add 53 I Iodine 4.93 g/cc
add 54 Xe Xenon 5.9 g/L
add 55 Cs Cesium 1.873 g/cc
add 56 Ba Barium 3.59 g/cc
add 57 La Lanthanum 6.15 g/cc
add 58 Ce Cerium 6.77 g/cc
add 59 Pr Praseodymium 6.77 g/cc
add 60 Nd Neodymium 7.01 g/cc
add 61 Pm Promethium 7.3 g/cc
add 62 Sm Samarium 7.52 g/cc
add 63 Eu Europium 5.24 g/cc
add 64 Gd Gadolinium 7.895 g/cc
add 65 Tb Terbium 8.23 g/cc
add 66 Dy Dysprosium 8.55 g/cc
add 67 Ho Holmium 8.8 g/cc
add 68 Er Erbium 9.07 g/cc
add 69 Tm Thulium 9.32 g/cc
add 70 Yb Ytterbium 6.9 g/cc
add 71 Lu Lutetium 9.84 g/cc
add 72 Hf Hafnium 13.31 g/cc
add 73 Ta Tantalum 16.65 g/cc
add 74 W Tungsten 19.35 g/cc
add 75 Re Rhenium 21.04 g/cc
add 76 Os Osmium 22.6 g/cc
add 77 Ir Iridium 22.4 g/cc
add 78 Pt Platinum 21.45 g/cc
add 79 Au Gold 19.32 g/cc
add 80 Hg Mercury 13.546 g/cc
add 81 Tl Thallium 11.85 g/cc
add 82 Pb Lead 11.35 g/cc
add 83 Bi Bismuth 9.75 g/cc
add 84 Po Polonium 9.3 g/cc
add 86 Rn Radon 9.73 g/L
add 88 Ra Radium 5.5 g/cc
add 89 Ac Actinium 10.07 g/cc
add 90 Th Thorium 11.724 g/cc
add 91 Pa Protactinium 15.4 g/cc
add 92 U Uranium 18.95 g/cc
add 93 Np Neptunium 20.2 g/cc
add 94 Pu Plutonium 19.84 g/cc
add 95 Am Americium 13.67 g/cc
add 96 Cm Curium 13.5 g/cc
add 97 Bk Berkelium 14.78 g/cc
add 98 Cf Californium 15.1 g/cc

all = add('get');
