 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = G * x	or x = G' * y
% Copyright 2002-2-20, Jeff Fessler, The University of Michigan

%
% scalar * G
%
if isa(ob, 'double') & length(x) == 1 & isa(x, 'Gtomo2_dsc')
	y = x;
	y.scale = ob;
	return
end

if ob.apower ~= 1, error notdone, end

%
% full forward projection
%
if ~ob.is_transpose

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

	y = wtfmex('dsc,proj', ob.arg', single(x), uint8(ob.mask), ...
		int32(0), int32(1), int32(ob.nthread), int32(ob.chat));


%
% full back-projection
%
else
	y = wtfmex('dsc,back', ob.arg', single(x), uint8(ob.mask), ...
		int32(0), int32(1), int32(ob.nthread), int32(ob.chat));

	if ob.is_masked
		if size(y,1) == numel(ob.mask)
			y = y(ob.mask,:);	% [nxy,nz] -> [np,nz]
		else
			y = y(ob.mask);		% [nx,ny] -> [np]
		end
	end
end

y = ob.scale * double(y);	% trick: no '(:)' here; let wtfmex do the job!
