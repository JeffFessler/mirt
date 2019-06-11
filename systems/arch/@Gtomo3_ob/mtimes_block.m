 function y = mtimes_block(ob, x, iblock, nblock)
%function y = mtimes_block(ob, x, iblock, nblock)
% y = G(i'th block) * x	or y = G'(i'th block) * x
% in either case the projection data will be "small"
% iblock is 1,...,nblock

% support 'exists' option for seeing if this routine is available
if nargin == 2 & ischar(x) & streq(x, 'exists')
	y = 1;
	return
end

if nargin ~= 4
	error(mfilename)
end

persistent Warned
if isempty(Warned), Warned = 0; end
if ~Warned
	warning 'slow: using full projections as a kludge for ordered subsets!'
	Warned = 1;
end

n1 = ob.n1;
n2 = ob.n2;
n3 = ob.n3;

%
% forward projection
%
if ~ob.is_transpose
	if ob.apower ~= 1, error notdone, end

	% if needed, expand concise column
	if ob.is_masked
		idim = size(x);
		np = sum(ob.mask(:));
		nxyz = numel(ob.mask);
		if idim(1) ~= np
			error 'size mismatch'
		end
		x = embed(x, ob.mask);
		x = reshape(x, nxyz, idim(2));
	end

        y = double(f3d_mex('proj', single(x), int32(ob.chat)));

	error 'not implemented!'

%	y = double(wtfmex('chat', ob.chat, 'proj,block', single(x), ...
%		int32(iblock-1), int32(nblock)));

if 0	% not implemented

	% fix: extract the relevent columns - should do in wtfmex!
	ia = iblock:nblock:na;

	if ndims(y) == 3	% [n1,n2,n3]
		y = y(:,ia,:);

	elseif size(y,1) == nb*na	% [nba,nz]
		nz = size(y,2);
		y = reshape(y, nb, na, nz);
		y = y(:,ia,:);
		y = reshape(y, nb*length(ia), nz);	% [n,nz]

	elseif size(y,1) == nb		% [nb,na]
		y = reshape(y, nb, na);
		y = y(:,ia);		% [nb,n]
	else
		error size
	end
end


%
% backprojection
%
else
	if ob.apower ~= 1, error notdone, end

error notdone

if 0
	ia = iblock:nblock:na;
	nv = length(ia);

	if ndims(x) == 3		% [n1,nv,nz]
		error todo
	elseif size(x,1) == n1		% [nb,nv]
		error todo
	elseif size(x,1) == nb*nv	% [nv,nz]
		nz = size(x,2);
		tmp = zeros(nb, na, nz);
		tmp(:,ia,:) = reshape(x, [nb nv nz]);
		x = reshape(tmp, [nb*na nz]);
	else
		error bug
	end

	y = double(wtfmex('chat', ob.chat, 'back,block', single(x), ...
		int32(iblock-1), int32(nblock)));

	if ob.is_masked
		if size(y,1) == numel(ob.mask)
			y = y(ob.mask,:);	% [nxy,nz] -> [np,nz]
		else
			y = y(ob.mask);		% [nx,ny] -> [np]
		end
	end
end

end
