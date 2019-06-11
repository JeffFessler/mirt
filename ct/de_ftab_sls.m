 function sls = de_ftab_sls(varargin)
%function sls = de_ftab_sls(varargin)
%|
%| Determine structure that characterizes the "s" limits for polyenergetic CT,
%| where s_l is a line integral through the lth material type, l=1...L.
%| L is the number of material components.
%|
%| in
%| option
%|	'sl'	cell{LL} sample thicknesses for each of LL materials	
%|	'n'	[L,1]	number of samples of lth material integrals
%|			default: [50 30] or [length(sl{ll})]
%|	'max'	[L,1]	maximum material density line integrals, units: g/cm^2
%|			default: [45 43] or [max(sl{ll})]
%|	'min'	[L,1]	minimum material density line integrals, units: g/cm^2
%|			default: [0] or [min(sl{ll})]
%|
%| out
%|	sls	strum
%|		data:
%|		sls.sl,n,max,min	see above
%|		methods:
%|		.sll			[n1,n2,...,nL,L] ndgrid of sl{*} values
%|
%| Copyright 2008-6-15, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test')
	de_ftab_sls_test
return
end
if ~nargout, help(mfilename), error(mfilename), end

sls.sl = {};
sls.n = [];
sls.min = [];
sls.max = [];

sls = vararg_pair(sls, varargin);

[sls.sl sls.n sls.min sls.max] = ...
	de_ftab_sls_do(sls.sl, sls.n, sls.min, sls.max);

sls = strum(sls, {'sll', @de_ftab_sls_sll, '()'});


% de_ftab_sls_sll()
function out = de_ftab_sls_sll(sls, varargin)
sll = ndgrid_jf('mat', sls.sl);
out = sll(varargin{:});


% de_ftab_sls_do()
function [sl, s_n, s_min, s_max] = de_ftab_sls_do(sl, s_n, s_min, s_max)

% number of samples of the material integrals "s"
if isempty(s_n)
	if isempty(sl)
		s_n = [45 43];
%		s_n = [51 31]; % makes plot_jac look nice
	else
		for ll=1:length(sl)
			s_n(1,ll) = length(sl{ll});
		end
	end
end

% minimum material "integrals"
if isempty(s_min)
	if isempty(sl)
		s_min = [0 0];
	else
		for ll=1:length(sl)
			s_min(1,ll) = min(sl{ll});
		end
	end
end

% maximum material "integrals"
if isempty(s_max)
	if isempty(sl)
		% soft max: 50cm * 1g/cc
		% bone max: 15cm * 2g/cc (for now)
		s_max = [50 30];
	else
		for ll=1:length(sl)
			s_max(1,ll) = max(sl{ll});
		end
	end
end

if isempty(sl)
	for ll=1:length(s_n)
		sl{ll} = linspace(s_min(ll), s_max(ll), s_n(ll))';
	end
end


% de_ftab_sls_test
function de_ftab_sls_test
sls = de_ftab_sls;
pr 'size(sls.sll)'
pr 'size(sls.sll(:,:,1))'
