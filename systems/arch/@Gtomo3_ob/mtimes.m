 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = G * x	or x = G' * y

%
% full forward projection
%
if ~ob.is_transpose
	if ob.apower ~= 1, error notdone, end

	% if needed, convert concise column to 3d array
	if ob.is_masked
		if size(x,1) ~= sum(ob.mask(:));
			error 'size mismatch'
		end
		x = embed(x, ob.mask);
	end

	x = single(reshape(x, [ob.nx ob.ny ob.nz]));
	y = double(f3d_mex('proj', x, int32(ob.chat)));


%
% full backprojection
%
else
	if ob.apower ~= 1, error notdone, end
	x = single(reshape(x, [ob.n1 ob.n2 ob.n3]));
	y = double(f3d_mex('back', x, int32(ob.chat)));
	if ob.is_masked
		y = y(ob.mask);
	end
end

y = y(:);
