 function y = mtimes_block(ob, x, iblock, nblock)
%function y = mtimes_block(ob, x, iblock, nblock)
% y = G(i'th block) * x	or y = G'(i'th block) * x
% in either case the project data will be "small"
% iblock is 1,...,nblock

% support 'exists' option for seeing if this routine is available
if nargin == 2 & ischar(x) & streq(x, 'exists')
	y = 1;
	return
end

if nargin ~= 4
	error(mfilename)
end

x = double(x); % sparse insists

if ob.apower == 1
	if ob.is_transpose
		y = (x' * ob.blocks{iblock})';
	else
		y = ob.blocks{iblock} * x;
	end

else

	B = ob.blocks{iblock};	% the ith block
	B = B .^ ob.apower;
	if ob.is_transpose
		B = B';
	end

	y = B * x;
end
