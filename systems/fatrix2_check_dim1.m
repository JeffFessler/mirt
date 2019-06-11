  function fatrix2_check_dim1(x, dim)
%|function fatrix2_check_dim1(x, dim)
%|
%| check that array x has dimensions "dim"
%| account for the "trailing 1" of size() for 1d column vectors

if nargin < 2, help(mfilename), error(mfilename), end

siz = size(x);
if numel(siz) == numel(dim) + 1 && siz(end) == 1
	siz = siz(1:end-1);
end
if ~isequal(siz, dim)
	fail('siz "%s" vs dim "%s" mismatch', num2str(size(x)), num2str(dim))
end
