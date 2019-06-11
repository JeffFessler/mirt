 function y = mtimes_block(ob, x, istart, nblock)
%function y = mtimes_block(ob, x, istart, nblock)
% y = G(i'th block) * x	or y = G'(i'th block) * x
% in either case the projection data will be "small"
% istart is 1,...,nblock

% support 'exists' option for seeing if this routine is available
if nargin == 2 && ischar(x) && streq(x, 'exists')
	y = ~isempty(ob.handle_mtimes_block);
	return
end

if nargin ~= 4
	error(mfilename)
end

if isempty(ob.handle_mtimes_block)
	error 'bug: no mtimes_block() method for this object'
end

if ob.is_transpose
	y = Fatrix_mtimes_block_back(ob, x, istart, nblock);
else
	y = Fatrix_mtimes_block_forw(ob, x, istart, nblock);
end

% caution: cascade_* may not work except for scalars here

function y = Fatrix_mtimes_block_forw(ob, x, istart, nblock);

x = do_cascade(ob.cascade_before, x, false, istart, nblock, true);
y = ob.handle_mtimes_block(ob.arg, ob.is_transpose, x, istart, nblock);
y = do_cascade(ob.cascade_after, y, false, istart, nblock, false);


function x = Fatrix_mtimes_block_back(ob, y, istart, nblock);

y = do_cascade(ob.cascade_after, y, true, istart, nblock, false);
x = ob.handle_mtimes_block(ob.arg, ob.is_transpose, y, istart, nblock);
x = do_cascade(ob.cascade_before, x, true, istart, nblock, true);
