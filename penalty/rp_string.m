 function [pstring, wstring, cstring] = rp_string(type, beta, delta)
%function [pstring, wstring, cstring] = rp_string(type, beta, delta)
%
% Specify roughness penalty string 
% pstring		for penalty
% wstring		for calculating penalty 'weights' (ala halfquad)
% cstring		for curvature (2nd derivative) of penalty
% beta and delta should be numeric values (full precision!)
%
% Copyright Mar. 1999, Jeff Fessler, The University of Michigan

if nargin < 2
	help rp_string
end

%
% Quadratic penalty (weights all unity)
% wstring = '%s = beta * ones(size(%s));'
%
if strcmp(type, 'quad')
	pstring = sprintf('%%s = %.18e * (%%s).^2 / 2;', beta);
	wstring = sprintf('%%s = %.18e * ones(size(%%s));', beta);
	cstring = wstring;

%
% Lange3 penalty
% wstring = '%s = beta ./ (1 + abs(%s) / delta);'
%
elseif strcmp(type, 'lange3')
	pstring = sprintf('%%s = lange3(%%s, %19.18e, %19.18e);', beta, delta);
	wstring = sprintf('%%s = %19.18e ./ (1 + abs(%%s) / %19.18e);', beta, delta);
	cstring = sprintf('%%s = %19.18e ./ (1 + abs(%%s) / %19.18e).^2;', beta, delta);

else
	error(['Unknown type "' type '"'])
end
