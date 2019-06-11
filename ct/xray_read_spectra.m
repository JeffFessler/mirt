  function xrs = xray_read_spectra(stype, varargin)
%|function xrs = xray_read_spectra(stype, [options])
%|
%| Read X-ray spectra data and initialize a structure that describes "M"
%| piecewise constant polyenergetic X-ray spectra, where M=2 for dual-kVp case.
%|
%| in
%|	stype		char	which spectrum model:
%|				'mono,60,100' dual mono-energetic
%|				'poly1,kvp1[,kvp2][,...]' polyenergetic
%|		or:	cell	{en [ne 1], sp [ne M]} sampled spectra
%|
%| option
%|	'en'		[ne 1]	specify energy sampling
%|	'filters'	cell{M} optional additional filtration of each spectrum
%|		{{mtype11, thick11, mtype12, thick12, ...}, ...
%|		 {mtypeM1, thickM1, mtypeM2, thickM2, ...}}
%|	'units'		''	'cm' (default) or 'mm' for filter thicknesses
%|	'show'		1|0	plot?
%|	'name'		{char}	names for each spectrum
%|	'kvp'		[M 1]
%|	'dir_spectra'	char	directory with spectra.  default:
%|				[path_find_dir('ct') filesep 'xray-spectra'];
%| out
%|	xrs		strum
%| data:
%|	xrs.en		[ne 1]	energy list
%|	xrs.sp		[ne M]	M spectra, for M kVp settings: I_m(E)
%|	xrs.Ide		[ne M]	differential array for "integrating": I_m(E) dE
%|	xrs.I		[M 1]	total intensity: sum(Ide)
%|	xrs.kvp		[M 1]	tube potentials (voltage in kVp)
%|	xrs.eff		[M 1]	mean energies (incident)
%| methods:
%|		.plot	plot spectra
%|
%| Note: xrs.eff perhaps should be named xrs.mean because it is a mean energy
%| of the incident spectrum.  It is not the "effective energy" as defined
%| e.g., on p56 of 4th edition of W. R. Hendee's "Medical Imaging Physics"
%| because that definition depends on the attenuation of a particular material.
%|
%| Copyright 2001-04-27, Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), return, end
if streq(stype, 'test'), xray_read_spectra_test, return, end

% defaults
arg.en = [];
arg.filters = {};
arg.show = false;
arg.name = {};
arg.kvp = [];
arg.units = 'cm';
arg.dir_spectra = [path_find_dir('ct') filesep 'xray-spectra'];
arg = vararg_pair(arg, varargin);

xrs.kvp = arg.kvp;

if ischar(stype)
	xrs = xray_read_spectra_char(stype, arg.en, arg.filters, arg.dir_spectra);
	xrs.type = stype;

elseif iscell(stype) && numel(stype) == 2
	xrs.en = stype{1};
	xrs.sp = stype{2};
	xrs.type = '';
	xrs.filters = arg.filters;
	if ~isequal(size(xrs.en,1), size(xrs.sp,1))
		error 'bad stype as cell usage'
	end
else
	error 'bad stype'
end

MM = size(xrs.sp, 2);

%
% apply filtration, if any
%
if ~isempty(xrs.filters)
	if numel(xrs.filters) ~= 1 && numel(xrs.filters) ~= MM
		error 'should be 1 or M sets of filters'
	end
	for mm=1:MM
		xrs.sp(:,mm) = xray_apply_filters(xrs.sp(:,mm), xrs.en, ...
			xrs.filters{min(mm,end)}, arg.units);
	end
end

%
% precompute some aspects of the spectra
%
xrs.Ide = xrs.sp .* repmat(difff(xrs.en), 1, MM); % [N M] I(E) dE
xrs.I = sum(xrs.Ide); % [1 M] spectrum integral
xrs.eff = xrs.en' * xrs.Ide ./ xrs.I;
for mm=1:MM
	xrs.sp_at_eff(mm) = interp1(xrs.en, xrs.sp(:,mm), xrs.eff(mm));
end

if isempty(arg.name)
	if ischar(stype)
		xrs.name = cat(2, {stype}, cell(1,MM-1));
	else
		xrs.name = '';
	end
else
	xrs.name = arg.name;
end
xrs = strum(xrs, {'plot', @xray_plot_spectra});

if arg.show
	xrs.plot;
end


%
% xray_read_spectra_char()
%
function xrs = xray_read_spectra_char(stype, en, filters, dir_spectra)
xrs.en = en;
xrs.filters = filters;

%
% simplest option is mono-energetic (for testing)
% usage: mono,kvp1,kvp2,...
%
if streq(stype, 'mono', 4)
	xrs.kvp = str2num(strrep(stype(6:end), ',', ' '));
	for mm=1:numel(xrs.kvp)
		xrs.en = [20:140]';
		xrs.sp(:,mm) = xrs.en == xrs.kvp(mm);
	end

%
% polyenergetic spectra
% usage: poly1,kvp1,kvp2,...
%
elseif streq(stype, 'poly1', 5)
	if ~exist(dir_spectra, 'dir')
		fail('spectra dir "%s" not in path', dir_spectra)
	end
	xrs.kvp = str2num(strrep(stype(7:end), ',', ' ')); % [M]

	% read raw data
	MM = numel(xrs.kvp);
	for mm=1:MM
		kvp = xrs.kvp(mm);
		raw = sprintf('spectra.%d', kvp);
		raw = [dir_spectra filesep raw];
		tmp = load_ascii_skip_header(raw);
		sr.enc{mm} = tmp(:,1);
		sr.spc{mm} = tmp(:,2);

		% The Wilderman/Sukovic spectra must be scaled by energy!
		sr.spc{mm} = sr.spc{mm} .* sr.enc{mm};
	end

	% interpolate onto same energy sampling
	[xrs.en xrs.sp] = xray_sp_interp(xrs.en, sr.enc, sr.spc);

%
% 1st spectra from predrag sukovic
%
elseif streq(stype, 'ps1')
	if ~exist(dir_spectra, 'dir')
		fail('spectra dir "%s" not in path', dir_spectra)
	end

	xrs.kvp = [80 140];

	for mm=1:numel(xrs.kvp)
		file = sprintf('xray%03d.mat', xrs.kvp(mm));
		file = [dir_spectra filesep file];
		if ~exist(file, 'file')
			error(sprintf('file "%s" not found', file))
		end
		raw = load(file);
		ie = raw.energy >= 20 & raw.energy <= 140;
		sr.enc{mm} = raw.energy(ie);
		sr.spc{mm} = raw.spe(ie) .* raw.energy(ie);
	end

	[xrs.en xrs.sp] = xray_sp_interp(xrs.en, sr.enc, sr.spc);

%
% spectra used for 2002 SPIE talk (fix: or did i just use 'ps1' ???)
%
elseif streq(stype, 'spie02')
	xrs = xray_read_spectra('poly1,80,140');
%	xrs.filters = { {{'aluminum', 0.25}, {'copper', 0.05}} };
	xrs.filters = {{'aluminum', 0.25, 'copper', 0.05}};
%	xrs.filters = {{'aluminum', 0.10}, {'copper', 0.05}};
%	xrs.filters = {{'aluminum', 0.5, 'copper', 0.04, 'csi', 0.1}};
%	xrs.filters = {{'aluminum', 0.5, 'copper', 0.04, 'gadolinium', 0.1}};


else
	fail('bad stype "%s"', stype)
end


%
% interpolate onto same energy sampling
%
function [en, sp] = xray_sp_interp(en, enc, spc)
MM = numel(spc);
if isempty(en)
	tmp = zeros(MM,1);
	for mm=1:MM
		tmp(mm) = max(enc{mm});
	end
	mm = imax(tmp); % find col with largest max
	en = enc{mm};
end
sp = zeros(numel(en), MM);
for mm=1:MM
	sp(:,mm) = interp1(enc{mm}, spc{mm}, ...
		en, 'linear', 0); % extrapolate with zeros
end


%
% xray_apply_filters()
%
function sp = xray_apply_filters(sp, en, filters, units);
% {mtype1, thick1, mtype2, thick2, ...}, ...
if ~iscell(filters)
	error 'filtration arguments must be cells'
end
for ii=1:2:numel(filters)
	sp = sp .* xray_filters(filters{ii}, filters{ii+1}, en, 'units', units);
end


%
% xray_plot_spectra()
% plot routine
%
function dummy = xray_plot_spectra(xrs, varargin)
dummy = [];

arg.title = xrs.type;
arg.kmin = min(xrs.en);
arg.kmax = max(xrs.en);
arg.spike = false;
arg.threshold = 1.2;
arg.hold = false;
arg.linetype = 'y.-';
arg = vararg_pair(arg, varargin);

MM = size(xrs.sp,2);
if ~im, return, end

if ~arg.hold, clf, end
pl=@(m) subplot(MM*100+10+m);
for mm=1:MM
	pl(mm)
	if arg.hold, hold on, end
	plot(xrs.en, xrs.sp(:,mm), arg.linetype, ...
		xrs.eff(mm) * [1 1], [0 xrs.sp_at_eff(mm)], 'm--')
	if arg.hold
		hold off
	else
		axis([arg.kmin arg.kmax [-0.00 1.05]*max(xrs.sp(:,mm))])
	end
	xt = [arg.kmin round(xrs.eff(mm)) arg.kmax];
	if arg.spike
		where = xrs_find_spike(xrs.sp(:,mm), arg.threshold);
		xt = [xt col(xrs.en(where))'];
	end
	xtick(sort(xt)), ytick([0])
	ylabelf('I%d(E)', mm)
	if mm == 1
		title(arg.title)
	end
	if ~isempty(xrs.name) && ~isempty(xrs.name{mm})
		title(xrs.name{mm})
	end

	if ~isempty(xrs.kvp)
		t = sprintf('%g kVp Spectrum', xrs.kvp(mm));
		texts(0.5, 0.8, t, 'color', 'green')
	end
end
xlabel 'Energy [keV]'


%
% xrs_find_spike()
% find spikes in spectrum
%
function where = xrs_find_spike(sp, threshold)
local = (sp(1:end-2) + sp(3:end)) / 2;
where = (sp(2:end-1) - local) > threshold * local;
where = 1 + find(where);


%
% xray_read_spectra_test()
% test routine
%
function xray_read_spectra_test

stype = 'mono,60,100';
stype = 'poly1,80,100';
stype = 'spie02';
stype = 'ps1';
stype = 'poly1,160,100,60';
filts = {{}, {}, {}};
xrs1 = xray_read_spectra(stype, 'filters', filts);
filts = {{'aluminum', 0.1}, {'copper', 0.2}, {}};
xrs2 = xray_read_spectra(stype, 'filters', filts);
xrs1.plot;
xrs2.plot('hold', true, 'linetype', 'c-');
%clf, plot(xrs1.en, xrs1.sp, '-', xrs2.en, xrs2.sp, '--')
