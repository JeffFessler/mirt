 function mas = xray_read_mac(mtype, varargin)
%function mas = xray_read_mac(mtype, [options])
%|
%| Read mass attenuation coefficients for a given material type.
%| Returns a strum that can interpolate onto desired energies.
%| 
%| in
%|	mtype		'aluminum', 'copper', 2, '2', '02-helium', ...
%|			See xray_material_file_name.m
%|			(Optionally can be a cell array of several materials.)
%|
%| options
%|	'units'		cm | mm	default: cm
%|	'interp'	{}	interpolator type.  default {'pchip', 'extrap'}
%|	'shortfile'	0|1	return short file name instead of full path
%|
%| out
%|	mas	strum
%| data:
%|	mas.type	{L}	material type
%|	mas.file	{L}	file for each material
%|	mas.units	'char'	cm | mm
%| methods:
%|	mas.plot		plot it
%|	mas.mac(kev)	[N L]	mass attenuation coefficients [cm^2/g],
%|				given vector of energies "kev"
%|	max.mean(kev, Ide) [M L] spectrally-weighted "mean" mac
%| private:
%|	mas.mac_raw	{L}	raw mass attenuation coefficients [cm^2/g]
%|	mas.mac_kev	{L}	raw energies from tables
%|
%| Copyright 2008-6-15, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(mtype, 'test'), xray_read_mac_test, return, end
if nargin == 1 && streq(mtype, 'test-water'), xray_read_mac_water, return, end

if isnumeric(mtype) && length(mtype) > 1
	mtype = num2cell(mtype);
end

if ~iscell(mtype)
	mtype = {mtype};
end

% defaults
mas.units = 'cm';
%arg.interp = {'linear', 'extrap'};
%arg.interp = {'spline', 'extrap'};
mas.interp = {'pchip', 'extrap'};
mas.shortfile = false;

mas = vararg_pair(mas, varargin);
mas.type = mtype;

is_mm = 0;
if streq(mas.units, 'mm')
	is_mm = 1;
elseif ~streq(mas.units, 'cm')
	error 'bad units'
end

if streq(mtype{1}, 'given') % special case {'given', mtype, mac, kev}
	mas.type = mtype{2};
	for ll=1:length(mas.type)
		mas.mac_raw{ll} = mtype{3}(:,ll);
		kev = mtype{4};
		if ncol(kev) > 1
			mas.kev_raw{ll} = kev(:,ll);
		else
			mas.kev_raw{ll} = kev(:,1);
		end
		clear kev
	end

elseif streq(mtype{1}, 'pca_do')
	fail 'not done'
%mac.type = {'pca1', 'pca2'}; % hardwired to first two components
%tmp = type{2}; % arguments for de_component2()
%tmp = de_component2('kev', en, tmp{:});
%mac.mac = tmp.basis; % [ne,2]

else
	for ll=1:length(mas.type)
		[mas.mac_raw{ll} mas.kev_raw{ll} mas.file{ll}] = ...
				xray_read_mac_raw(mas.type{ll}, is_mm);

		if mas.shortfile
			mas.file{ll} = regexprep(mas.file{ll}, '.*\/', '');
		end
	end
end

meth = {'plot', @xray_read_mac_plot, '()';
	'mac', @xray_read_mac_interp, '([kev])';
	'mean', @xray_read_mac_mean, '(kev, Ide)';
	};
mas = strum(mas, meth);


%
% xray_read_mac_raw()
%
function [mac kev file] = xray_read_mac_raw(mtype, is_mm)

file = xray_material_file_name(mtype);
tmp = load_ascii_skip_header(file); % read uncommented lines
kev = tmp(:,1) * 1000;	% keV
mac = tmp(:,2);	% mass atten. coeff. [cm^2/g]

if is_mm
	mac = mac * 100;
end

if 0 % look at k-edges
	jump = find(diff(kev) == 0);
	if ~isempty(jump)
		file = regexprep(file, '.*\/', '');
		jump = num2str(kev(jump)', '%5.1f');
		printf('%s %s', file, jump)
	end
end


%
% xray_read_mac_interp()
% interpolate onto desired energies
%
function mac = xray_read_mac_interp(mas, kev)

if ~isvar('kev') || isempty(kev)
	kev = mas.kev{1};
	if length(mas.type) ~= 1
		warn('using energies of 1st material type')
	end
end

LL = length(mas.type);
mac = zeros(length(kev), LL); % [N,L]

% trick: allow for the k-edge jumps!
interp = mas.interp;
for ll = 1:LL
	tmp = log(mas.mac_raw{ll}); % interpolate on a log scale
	mac(:,ll) = interp1_jump(mas.kev_raw{ll}, tmp, kev, ...
		'monodown', interp{:});
end
mac = exp(mac);


%
% xray_read_mac_mean()
% spectrally-weighted "mean" mac
%
function bar = xray_read_mac_mean(mas, kev, Ide)
mac = mas.mac(kev);
bar = diag(1 ./ sum(Ide)) * (Ide' * mac);


%
% xray_read_mac_plot()
%
function dummy = xray_read_mac_plot(mas, varargin)
dummy = [];

arg.kev = [10:510]';
arg = vararg_pair(arg, varargin);
kev = arg.kev;

mac = mas.mac(kev);

argp = {};
for ll=1:length(mas.type)
	kev_raw = mas.kev_raw{ll};
	ie = min(kev) <= kev_raw & kev_raw <= max(kev);
	argp = {argp{:}, kev_raw(ie), mas.mac_raw{ll}(ie), 'o'};
end

semilogy(argp{:})
mtypes = mas.type;
legend(mtypes{:})
hold on, semilogy(kev, mas.mac(kev), '-'), hold off
xlabel 'KeV', ylabel 'mass atten. coef.', axis tight


%
% xray_read_mac_test()
% example usage
%
function xray_read_mac_test
mtypes = {'lead', 'aluminum', 'water', 'lexan'};
mtypes = {'water', 'bone', 'iodine', 'adipose'};
%mtypes = {'water', 'bone', 'calcium'};
%mtypes = {'04-beryllium', 'water'}
%mtypes = {'iodine', 'cesium', 'csi'};
%mtypes = {'cadmium', 'tellurium', 'cdte'};
%mtypes = {'06', '06-carbon-graphite', 'carbon-graphite'};
mas = xray_read_mac(mtypes);
mas.plot('kev', 20:10:2000);

xrs = xray_read_spectra('ps1');
mas.mean(xrs.en, xrs.Ide);

if 0 % examine CdTe
	% Cd atomic mass about 112
	% Te atomic mass about 128
	en = 5:500;
	tmp = mas.mac(en);
	tmp = [tmp(:,1) * 112 / (112 + 128) + tmp(:,2) * 128 / (112 + 128) ...
		tmp(:,3)];
	plot(en, tmp(:,1) - tmp(:,2))
end


% xray_read_mac_water()
% test mac for water
%
function xray_read_mac_water
wmas = xray_read_mac({'water'});
emas = xray_read_mac({'hydrogen', 'oxygen'});
kev = 40:1:200;
%emas.plot('kev', kev);
tmp = emas.mac(kev);
pred = 2/18 * tmp(:,1) + 16/18 * tmp(:,2); % approximate mac
plot(kev, pred, '-', kev, wmas.mac(kev), '--')
legend('predict', 'nist')
%plot(kev, pred-wmas.mac(kev), '--')
