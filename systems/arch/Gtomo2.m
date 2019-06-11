 function ob = Gtomo2(chat)
%function ob = Gtomo2(chat)
%
%	create default "parent" 2D tomography object (structure)
%	with some basic methods common to all child tomography objects
%	chat: 0 to disable debugging messages, >=1 to display them
%
%	Copyright Mar 2001, Jeff Fessler, University of Michigan

if nargin < 1, chat = 0; end

%
%	the dimensions are usually [ob.nb * ob.na, ob.nx * ob.ny]
%	unless masked, in which case [ob.nb * ob.na, sum(mask(:))]
%	or if tranposed (in which case swapped)
%
ob.dims = [0 0];

ob.nx = 0;	% input image size
ob.ny = 0;
ob.nb = 0;	% output sinogram size (bins)
ob.na = 0;	% output sinogram size (angles)

ob.apower	= 1;		% array power, becomes 2 for G.^2
ob.scale	= 1;		% global scaling factor used in some children
ob.mask		= [];		% [nx,ny] logical array
ob.chat		= chat;
ob.index1	= [];		% row indices (unless transposed) e.g. G(3:9,:)
ob.index2	= [];		% col indices (unless transposed) e.g. G(:,3:9)
				% the empty default means *all* rows/cols
ob.is_transpose = false;
ob.is_masked	= false;	% set to 1 if G(:,mask(:))
%ob.version	= 1.0;

%ob = class(ob, 'Gtomo2');	% matlab inheritence sucks
