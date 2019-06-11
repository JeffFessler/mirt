  function [out ii] = jf_pair_parse(cells, str, varargin)
%|function [out ii] = jf_pair_parse(cells, str, varargin)

arg.suffix = ' ';
arg = vararg_pair(arg, varargin);

str = [str arg.suffix];

tf = strncmp(str, cells, numel(str));
ii = find(tf);
if numel(ii) < 1
	fail('no match for "%s"', str)
end
if numel(ii) > 1
	fail('%d matches for "%s"', numel(ii), str)
end

out = cells{ii};
out = out((numel(str)+1):end);
