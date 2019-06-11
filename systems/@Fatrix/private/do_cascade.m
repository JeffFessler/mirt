 function y = do_cascade(cascade, x, is_transpose, istart, nblock, is_before)
%function y = do_cascade(cascade, x, is_transpose, istart, nblock, is_before)
% do cascade * x or cascade' * x

if isempty(cascade)
	y = x;

elseif isa(cascade, 'function_handle') || isa(cascade, 'inline')
	if nargin(cascade) == 4
		y = cascade(x, is_transpose, istart, nblock);
	elseif nargin(cascade) == 2 || nblock == 1
		y = cascade(x, is_transpose);
	else
		error 'cascade_* should have 2 or 4 arguments?'
	end

else % matrix or object
	if is_transpose
		if is_before
			do_cascade_warn(cascade, nblock)
		end
		y = cascade' * x;
	else
		if ~is_before
			do_cascade_warn(cascade, nblock)
		end
		y = cascade * x;
	end
end

function do_cascade_warn(cascade, nblock)
if max(size(cascade)) > 1 && nblock ~= 1
	warning 'non scalar array/object cascade not tested with block!'
end
