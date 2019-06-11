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

nb = ob.nb;
na = ob.na;

%
% forward projection
%
if ~ob.is_transpose
	if ob.apower ~= 1, error notdone, end

	% if needed, expand concise column
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

	y = double(wtfmex('dsc,proj', ob.arg', single(x), uint8(ob.mask), ...
		int32(iblock-1), int32(nblock), ...
		int32(ob.nthread), int32(ob.chat)));

	% fix: extract the relevent columns - should do in wtfmex?
	ia = iblock:nblock:na;
	nv = length(ia);

	if ndims(y) == 3			% [nb,na,nz]
		y = y(:,ia,:);			% [nb,nv,nz]

	elseif size(y,1) == nb*na		% [nba,nz]
		nz = size(y,2);
		y = reshape(y, nb, na, nz);	% [nb,na,nz]
		y = y(:,ia,:);			% [nb,nv,nz]
		y = reshape(y, nb*nv, nz);	% [nb*nv,nz]

	elseif size(y,1) == nb			% [nb,na]
		y = y(:,ia);			% [nb,nv]

	else
		error size
	end


%
% backprojection
%
else
	if ob.apower ~= 1, error notdone, end

	ia = iblock:nblock:na;
	nv = length(ia);

	if ndims(x) == 3		% [nb,nv,nz]
		error todo

	elseif size(x,1) == nb*nv	% [nb*nv,nz]
		nz = size(x,2);
		tmp = zeros(nb, na, nz);
		tmp(:,ia,:) = reshape(x, [nb nv nz]);
		x = reshape(tmp, [nb*na nz]);

	elseif size(x,1) == nb		% [nb,nv]
		error todo

	else
		error bug
	end

	y = double(wtfmex('dsc,back', ob.arg', single(x), uint8(ob.mask), ...
		int32(iblock-1), int32(nblock), ...
		int32(ob.nthread), int32(ob.chat)));

	if ob.is_masked
		if size(y,1) == numel(ob.mask)
			y = y(ob.mask,:);	% [nxy,nz] -> [np,nz]
		else
			y = y(ob.mask);		% [nx,ny] -> [np]
		end
	end
end
