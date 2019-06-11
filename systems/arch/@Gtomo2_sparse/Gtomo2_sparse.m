 function ob = Gtomo2_sparse(file, chat, nx, ny, nb, na, mask)
%function ob = Gtomo2_sparse(file[, chat, nx, ny, nb, na, mask])
%
% Construct Gtomo2_sparse object, either from a .wtf file,
% or from a Matlab sparse matrix.
% This just uses an ordinary Matlab sparse matrix for the core!
% See Gtomo2_sparse_test.m for example usage.
% You create an system object by calling:
%	G = Gtomo2_sparse(file)
% and then you can use it thereafter by typing commands like
%	y = G * x.
% If the input is a sparse matrix, then give the dimensions nx,ny,nb,na.
% If the file argument is 'identity', then make sparse identity matrix case.
%
% If the input matrix is already "masked," then the mask argument is key.
%
% Copyright 2001-1-30, Jeff Fessler, The University of Michigan

if nargin < 2, chat = 0; end

%
% create default object, as required by Mathworks
%
ob = Gtomo2(chat);
ob.file = '';
ob.G = [];
ob.blocks = {};	% place to store blocks of G for subset algorithms

if nargin == 0	% required by Mathworks
%	help(mfilename)
%	warning 'Gtomo2_sparse called with no arguments!?'
	ob = class(ob, 'Gtomo2_sparse');
return
end


%
% trick for sparse identity matrices
%
if ischar(file) & streq(file, 'identity')
	if nb ~= nx | na ~= ny, error 'bug', end
	file = speye(nx*ny);
	file = file(:,mask(:)); % this is now a sparse matrix!
end


%
% if input is a .wtf file
%
if ischar(file)
	ob.file = file;

	[ob.G ob.nx ob.ny ob.nb ob.na] = wtfmex('load', file);
	ob.dims = [ob.nb * ob.na, ob.nx * ob.ny];

	tmp = sum(ob.G) > 0;
	ob.mask = reshape(tmp, ob.nx, ob.ny);

elseif issparse(file)
	ob.G = file;
	if nargin < 6, help(mfilename), error(mfilename), end
	ob.nx = nx;
	ob.ny = ny;
	ob.nb = nb;
	ob.na = na;
	ob.dims = [ob.nb * ob.na, ob.nx * ob.ny];
	if nb*na ~= size(ob.G,1)
		error 'bad row dimension'
	end
	if nx*ny ~= size(ob.G,2)
		if ~isvar('mask'), error 'need mask', end
		if any([nx ny] ~= size(mask)), error 'bad mask size', end
		ob.dims(2) = sum(mask(:));
	end
	if isvar('mask')
		ob.mask = mask;
	end
else
	error 'input must be filename or sparse matrix'
end

%ob.is.empty = false;
ob = class(ob, 'Gtomo2_sparse');
