 function [ob, args] = vararg_get(ob, key, varargin)
%function [ob, args] = vararg_get(ob, key, varargin)
% This function is used for parsing name/value pair arguments.
% For example, a function is called via fun( ..., 'key1', value1, ...) 
% and we want to add the value to the structure via: ob.key1 = value1;
if nargin < 1, ir_usage, end

if ~length(varargin)
	error(sprintf('need value for %s', key))
end

arg = varargin{1};
ob.(key) = arg;
%ob = setfield(ob, key, arg);
args = {varargin{2:end}};

if isfield(ob, 'chat') && ob.chat
	if isnumeric(arg)
		if max(size(arg)) == 1	% only print scalars
			printf('\t%s = %g', key, arg)
		else	% otherwise print sum (e.g., for mask)
			printf('\tsum(%s) = %g', key, sum(arg(:)))
		end

	elseif isstruct(arg)	% show first structure element, if string
		fields = fieldnames(arg);
%		name = getfield(arg, fields{1});
		name = arg.(fields{1});
		if ischar(name)
			printf('\t%s.%s = %s', key, fields{1}, name)
		end

	elseif ischar(arg)
		printf('\t%s = %s', key, arg)

	elseif isa(arg, 'cell')
		printf('\t%s {cell}', key)

	else
		warning 'display type not implemented'
	end
end
