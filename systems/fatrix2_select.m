  function z = fatrix2_select(mask, dim, y)
%|function z = fatrix2_select(mask, dim, y)
%|
%| version of masker that handles an empty mask
%|
%| in
%|	y	[*dim L]
%|	dim	[]		size(mask) if mask is nonempty
%|	mask	[(dim)]		possibly empty
%|
%| out
%|	z	[N L]		N = sum(mask(:))
%|
%| Copyright 2010-12-21, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

if ~isempty(mask)
	z = y(mask(:),:);
else
	z = y;
end
