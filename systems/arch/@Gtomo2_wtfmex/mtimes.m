 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = G * x	or x = G' * y
% Copyright 2002-2-20, Jeff Fessler, The University of Michigan

%
% scalar * G
%
if isa(ob, 'double') & length(x) == 1 & isa(x, 'Gtomo2_wtfmex')
	y = x;
	y.scale = ob;
	return
end

%
% full forward projection
%
if ~ob.is_transpose

	if ob.apower ~= 1, error notdone, end

	% if needed, expand concise column(s)
	if ob.is_masked
		idim = size(x);
		np = sum(ob.mask(:));
		nxy = numel(ob.mask);
		if idim(1) ~= np
			error 'size mismatch'
		end
		x = embed(x, ob.mask);
		x = reshape(x, nxy, idim(2));
	end

	y = wtfmex('chat', ob.chat, 'mult', single(x));


%
% full back-projection
%
else

	if ob.apower == 1
		y = wtfmex('chat', ob.chat, 'back', single(x));
	elseif ob.apower == 2
		y = wtfmex('chat', ob.chat, 'back2', single(x));
	else
		error 'only power 1,2 done'
	end

	if ob.is_masked
		if size(y,1) == numel(ob.mask)
			y = y(ob.mask,:);	% [nxy,nz] -> [np,nz]
		else
			y = y(ob.mask);		% [nx,ny] -> [np]
		end
	end
end

y = ob.scale * double(y);	% trick: no '(:)' here; let wtfmex do the job!
