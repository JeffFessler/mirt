 function iwater = de_ftab_iwater(s_arg, varargin)
%function iwater = de_ftab_iwater(s_arg, [options])
%
% Build "iwater" object that handles inverse of water-only beam hardening.

todo: under construction!

%
% in:	(these all have sensible defaults, so try calling with no args)
%	s_arg	cell	arguments for xray_read_spectra
% option
%	'mtype	cell	material types, e.g., {'soft', 'bone'}
%	'ftype	char	fit type for de_ftab_fit()
%	'sl	cell{LL} sample thicknesses for each of LL materials	
%	's_n	[L,1]	number of samples of l'th material integrals
%	's_max	[L,1]	maximum material density line integrals, units: g/cm^2
%	'wt_fm	{M}	weighting for fitting fm
%
% out:
%	ftab	struct
%	methods:
%		ftab.fit.fmfun(s1, s2, ...)
%
% Copyright 2001-04-27, Jeff Fessler, The University of Michigan
if nargin < 1, help(mfilename), error(mfilename), end

ftab.show = false;
ftab.sl = {};
ftab.s_n = [];
ftab.s_max = [];
ftab.wt_fm = {};
ftab.stop_after_fit = false;
ftab.mtype = {'soft', 'bone'};	% default materials
ftab.ftype = ''; % defer to de_ftab_fit() default
ftab = vararg_pair(ftab, varargin);

%
% Build inverse of soft-tissue component
% for conventional "water only" beam-hardening correction.
%
if ~isvar('ftab.inv_water'), printm 'water inv'

	fw = ftab.feval(ftab, ftab.sl{1}, 0);	% mm=1 spectrum, water only
	fw = squeeze(fw(:,1,:));	% [s_n(1) M]
	ftab.inv_water_f = linspace(0, max(fw(:)), 51)';
	for mm=1:MM
		ftab.inv_water(:,mm) = interp1(fw(:,mm), ftab.sl{1}, ...
			ftab.inv_water_f, 'cubic', 'extrap');
	end
	if any(isnan(ftab.inv_water)), error 'nan', end

	ftab.inv_water_eval = @(ftab, m, f) interp1(ftab.inv_water_f, ftab.inv_water(:,m), f, 'cubic', 'extrap');

	if ftab.show
		clf, subplot(221)
		plot(ftab.sl{1}, fw(:,1), '-', ...
			ftab.inv_water(:,1), ftab.inv_water_f, '.')
		axis tight, xlabel 's1', ylabel 'f1', title 'water fit m=1'
		subplot(222)
		plot(ftab.sl{1}, fw(:,2), '-', ...
			ftab.inv_water(:,2), ftab.inv_water_f, '.')
		axis tight, xlabel 's1', ylabel 'f2', title 'water fit m=2'

		wat(:,1) = ftab.inv_water_eval(ftab, 1, fw(:,1));
		wat(:,2) = ftab.inv_water_eval(ftab, 2, fw(:,2));
		err = wat - [ftab.sl{1} ftab.sl{1}];

		subplot(223)
		plot(ftab.sl{1}, err(:,1), '.-')
		axis tight, xlabel s1, ylabel err1, title 'Water inverse error'
		subplot(224)
		plot(ftab.sl{1}, err(:,2), '.-')
		axis tight, xlabel s1, ylabel err2, title 'Water inverse error'

		prompt
	end
	clear fw t2 wat err
end
