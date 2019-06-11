 function ob = Gnufft(varargin)
%function ob = Gnufft([mask,] args)
% Construct Gnufft object, which computes nonunform DSFT samples
% of signals with dimensions [(Nd)] via the NUFFT.
%
% The arguments are simply a cell array of the all the arguments
% that will be passed to "nufft_init()" in the appropriate order.
% See Gnufft_test.m for example usage.
%
% Basically, you create a system matrix object by calling:
%	G = Gnufft( ... )
% and then you can use it thereafter by typing commands like
%	y = G * x;
% which will auto-magically evaluate the DSFT samples.
% This is useful for iterative image reconstruction in MRI.
%
% Besides simple utilities like display, there are the following
% capabilities of this object:
%	y = G * x		forward operation
%	x = G' * y		adjoint operation
%
% Optional arguments
%	mask		logical support array
%
% Copyright 2003-6-1, Jeff Fessler, The University of Michigan

%
% default object
%
ob.dims = [0 0];
ob.apower	= 1;		% array power, becomes 2 for G.^2
%ob.scale	= 1;		% global scaling factor used in some children
ob.mask		= [];		% [(Nd)] logical array
ob.chat		= false;
ob.index1	= [];		% row indices (unless transposed) e.g. G(3:9,:)
ob.index2	= [];		% col indices (unless transposed) e.g. G(:,3:9)
				% the empty default means *all* rows/cols
ob.is_transpose = false;
ob.is_masked	= false;	% set to 1 if G(:,mask(:))
ob.is_subref	= false;	% set to 1 if subscripted

ob.arg = {};
ob.st = [];
ob.Nd = [];

if nargin < 1
	warning 'Gnufft called with too few arguments'
	help Gnufft
	ob = class(ob, 'Gnufft_ob');
	return

elseif nargin == 1
	% kludge to allow re-construction of object from structure
	if isstruct(varargin{1})
		ob = class(varargin{1}, 'Gnufft_ob');
		return
	end
	ob.mask = [];	% true(nx,ny);
	args = varargin{1};

elseif nargin == 2
	ob.mask = varargin{1};
	if ~islogical(ob.mask), error 'mask must be logical', end
	args = varargin{2};

else
	error 'too many arguments'
end

if ~isa(args, 'cell'), error 'args must be cell array', end

omega = args{1};
ob.Nd = args{2};
if length(ob.Nd) == 1, ob.Nd = [ob.Nd 1]; end	% handle 1D case
ob.dims = [size(omega,1) prod(ob.Nd)];

%
% initialize the structure needed by NUFFT
%
ob.arg = args;
ob.st = nufft_init(args{:});

% if mask is provided, then it is masked
if ~isempty(ob.mask)
	ob.is_masked = 1;
	ob.dims(2) = sum(ob.mask(:));
end

%ob.is.empty = false;
ob = class(ob, 'Gnufft_ob');
