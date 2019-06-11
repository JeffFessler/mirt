 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = G * x	or x = G' * y

%
% full forward projection
%
if ~ob.is_transpose

	% if needed, convert 2d array to concise column
	if ob.is_masked & size(x,1) ~= ob.dims(2);
		error 'not done'
%		if size(x,1) ~= sum(sum(ob.mask));
%			error 'size mismatch'
%		end
%		x = embed(x, ob.mask);
%		x = x(:);
	end

	if ob.apower ~= 1
		y = (ob.G).^(ob.apower) * x;
	else
		y = ob.G * x;
	end


%
% full backprojection
%
else
	if ob.apower ~= 1
		y = (ob.G .^ ob.apower)' * x;
	else
		y = (x' * ob.G)';
	end
%	if ob.is_masked
%		y = y(ob.mask,:);
%	end
end
