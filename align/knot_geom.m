 function st = knot_geom(varargin)
%function st = knot_geom(varargin)
%
% Create a "knot geometry" structure that describes the B-spline
% grid of a single 2d image depending on "image geometry"
% Using this structure should facilitate "object oriented" code.
%
% options for 2d
%	'nx'		# knots of x
%	'ny'		# knots of y (default: nx)
%	'mx'		scale in the x direction
%	'my'		scale in the y direction (default: mx)
%	'offset_x'	unitless (default: 0, leftmost)
%	'offset_y'	unitless (default: 0, uppermost)
%	'mask'		logical support mask
%
% options for 3d
%	'nz'		# knots of z
%	'mz'		scale in the z direction
%	'offset_z'	unitless (default: 0)
%
% out:
%	st	(strum)	initialized structure with methods
%
% methods:
%	st.x	x coordinates of each pixel
%	st.y	y coordinates of each pixel 
%	st.np	sum(st.mask(:)) (# of pixels to be estimated)
%	st.embed(column)	turn short column(s) into array(s)
%	st.maskit(x)		the opposite of embed
%	st.shape(x)		reshape long column x to nx,ny array
%
% methods for 3d:
%	st.z		z coordinates of each pixel
%
% Copyright October 2006, Se Young Chun and Jeff Fessler, University of Michigan
if nargin < 1, help(mfilename), error(mfilename), end

% defaults
st.nx = [];
st.ny = [];
st.nz = [];
st.mx = [];
st.my = [];
st.mz = [];
st.offset_x = 0;
st.offset_y = 0;
st.offset_z = 0;
st.mask = [];

st = vararg_pair(st, varargin);

% dimensions
if isempty(st.ny), st.ny = st.nx; end
if isempty(st.my), st.my = st.mx; end

if isempty(st.mz) 
	st.is3 = false;
	st.dim = [st.nx st.ny];
	st = rmfield(st, {'nz', 'mz', 'offset_z'});
else
	st.is3 = true;
	st.dim = [st.nx st.ny st.nz];
end

% mask
if isempty(st.mask)
	st.mask = true(st.dim);
elseif ~islogical(st.mask)
	error 'mask must be logical'
elseif ndims(st.mask) ~= 2 + st.is3 ...
	|| size(st.mask,1) ~= st.nx ...
	|| size(st.mask,2) ~= st.ny ...
	|| (st.is3 && size(st.mask,3) ~= st.nz)
	size(st.mask), st.nx, st.ny
	error 'bad input mask size'
end

meth = { ...
	'embed', @knot_geom_embed, '()'; ...
	'shape', @knot_geom_shape, '()'; ...
	'maskit', @knot_geom_maskit, '()'; ...
	'x', @knot_geom_x, '()'; ...
	'y', @knot_geom_y, '()'; ...
	'np', @knot_geom_np, '()'; ...
	};

if st.is3
	meth(end+[1],:) = { ...
		'z', @knot_geom_z, '()'; ...
		};
end
st = strum(st, meth);



% knot_geom_np()
function np = knot_geom_np(st)
np = sum(st.mask(:));



% knot_geom_embed()
function x = knot_geom_embed(st, x)
if issparse(x), x = knot_geom_embed_sparse(st, x); return, end
x = embed(x, st.mask);



% knot_geom_embed_sparse()
% trick: this is for "unpacking" sparse system matrices
function x = knot_geom_embed_sparse(st, x)
[i j a] = find(x);
ind = find(st.mask);
j = ind(j);
x = sparse(i, j, a, size(x,1), prod(st.dim));



% knot_geom_maskit()
% opposite of embed
% in 3d case, if input is [nx ny nz n4], output is [np n4] where np = sum(mask)
function x = knot_geom_maskit(st, x)
dim = size(x);
x = reshape(x, prod(st.dim), []);
x = x(st.mask(:),:);
if length(dim) > length(st.dim)
	x = reshape(x, [], dim((1+length(st.dim)):end));
end



% knot_geom_shape()
function x = knot_geom_shape(st, x)
if st.is3
	x = reshape(x, st.nx, st.ny, st.nz, []);
else
	x = reshape(x, st.nx, st.ny, []);
end



% knot_geom_x()
function x = knot_geom_x(st, varargin)
x = [st.offset_x : st.mx : st.mx*(st.nx-1) + st.offset_x];
x = x(varargin{:});



% knot_geom_y()
function y = knot_geom_y(st, varargin)
y = [st.offset_y : st.my : st.my*(st.ny-1) + st.offset_y];
y = y(varargin{:});



% knot_geom_z()
function z = knot_geom_z(st, varargin)
z = [st.offset_z : st.mz : st.mz*(st.nz-1) + st.offset_z];
z = z(varargin{:});

