 function y = mtimes_block(ob, x, iblock, nblock)
%function y = mtimes_block(ob, x, iblock, nblock)
%	y = G(i'th block) * x	or y = G'(i'th block) * x
%	in either case the project data will be "small"
%	iblock is 1,...,nblock

% support 'exists' option for seeing if this routine is available
if nargin == 2 & ischar(x) & streq(x, 'exists')
	y = 1;
	return
end

if nargin ~= 4
	error(mfilename)
end

if nblock ~= 1 | iblock ~= 1
	error 'no block for mtimes_block'
end
y = ob * x;
