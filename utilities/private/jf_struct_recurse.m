 function out = jf_struct_recurse(pn, st, varargin)
%function out = jf_struct_recurse(pn, st, varargin)
%|
%| Recursively descend a structure and examine its (numeric) members
%|
%| Copyright 2010-04-05, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

arg.type = 'minmax_nan';
arg.other = false; % show other types (e.g., char?)
arg.top = '';
arg.fun = [];
arg = vararg_pair(arg, varargin);

switch arg.type
case 'minmax_nan'
	arg.fun = @(x) minmax_nan(x);
otherwise
	if isempty(arg.fun)
		fail 'need ''type'' or ''fun'' option'
	end
end

jf_struct_recurse_do(st, arg.top, arg.fun, arg)
out = [];


% jf_struct_recurse_do()
function jf_struct_recurse_do(x, prefix, fun, arg)

%if isa(x, 'strum')
%	x = struct(x);
%end

if isempty(x)
	% ignore

elseif isnumeric(x) || islogical(x)
	r = fun(x);
	if numel(r) == 1
		printm('%s : %g', prefix, r)
	elseif numel(r) == 2
		printm('%s : %g %g', prefix, r(1), r(2))
	else
		printm('%s :', prefix)
		disp(r)
	end

elseif ischar(x)
	if arg.other
		printm('%s : %s', prefix, x)
	end

elseif iscell(x)
	for ii = 1:length(x)
		tmp = [prefix sprintf('{%d}', ii)];
		jf_struct_recurse_do(x{ii}, tmp, fun, arg)
	end

elseif isstruct(x)
	for ix = 1:numel(x) % in case x is a struct array
		y = x(ix);
		names = fieldnames(y);
		for ii = 1:length(names)
			tmp = [prefix '.' names{ii}];
			jf_struct_recurse_do(y.(names{ii}), tmp, fun, arg)
		end
	end

elseif isa(x, 'function_handle')
	% ignore

else
	try % try to convert other classes to struct, e.g., strum or fatrix
		x = struct(x);
		jf_struct_recurse_do(x, prefix, fun, arg)
	catch
		warn('unknown class %s', class(x))
	end
end


function out = min_nan(x)
if any(isnan(x(:)))
	out = nan;
else
	out = min(x(:));
end


function out = max_nan(x)
if any(isnan(x(:)))
	out = nan;
else
	out = max(x(:));
end


function out = minmax_nan(x)
out = [min_nan(x) max_nan(x)];
