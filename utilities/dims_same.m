  function y = dims_same(a, b, varargin)
%|function y = dims_same(a, b, varargin)
%| see if two arrays have the same dimensions
%| options
%|	'scalar_ok'	1 means "scalar b matches anything"
%|	'up_to_dim'	only match up to this dimension

if nargin < 2, ir_usage, end

arg.scalar_ok = 0;
arg.up_to_dim = [];
arg = vararg_pair(arg, varargin);
%if isempty(arg.up_to_dim)
%	arg.up_to_dim = max(ndims(a), ndims(b));
%end

if arg.scalar_ok && numel(b) == 1
	y = 1;
	return
end

if ~isempty(arg.up_to_dim)
	if ndims(a) < arg.up_to_dim || ndims(b) < arg.up_to_dim
		y = 0;
		return
	end
	asize = size(a);
	bsize = size(b);

	if any(asize(1:arg.up_to_dim) ~= bsize(1:arg.up_to_dim))
		y = 0;
		return
	end

else
	if ndims(a) ~= ndims(b)
		y = 0;
		return
	end

	if any(size(a) ~= size(b))
		y = 0;
		return
	end
end

y = 1;
