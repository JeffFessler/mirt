  function ob = Cdiff(kappa, varargin)
%|function C = Cdiff(kappa, [options])
%|
%| This object has been made obsolete by Cdiff1.m and Reg1.m
%|
%| Construct Cdiff object that can compute C * x and the adjoint C' * d
%| for a "differencing" matrix C for roughness penalty regularization.
%| This object is used internally by the roughness penalty object Robject,
%| which is designed for the form R(x) = \sum_k w_k \pot([Cx]_k)
%|
%| By construction, C = Diag(scale) * C1, where the elements of C1 are
%| all zeros except for a +1 and -1 in each row, and scale = sqrt(w_k)
%| where w_k = kappa_j * kappa_k / distance^distance_power.
%|
%| Caution: because the default 'edge_type' is 'none', this differencing
%| matrix will cause wrap-around effects at the image edges.  The user
%| must multiply this by a suitable diagonal matrix to avoid it.
%| This object is intended to be used "internally" by other objects
%| like Robject that provide the suitable diagonal matrix.
%|
%| in
%|	kappa	[nx ny nz]	binary support mask, or the "kappa" array,
%|				or just a vector giving dimensions like [64 64]
%|				If kappa is binary, then each row of C will
%|				be scaled by 1 or 0 or 1/sqrt(2) etc. 
%|				If kappa is floats, then each row will be
%|				scaled by sqrt(kappa_j kappa_k).
%| options
%|	'offsets' [1 nx ...]	manually specify neighbor offsets
%|				default in 2D is: [1 nx nx+1 nx-1];
%|				or a string option:
%|					'2d,hvd'	8 nearest neighbors
%|	'dosparse' [1 or 0]	return sparse matrix rather than Fatrix object
%|				(default is 0)
%|	'edge_type' [type]	edge conditions:
%|		'tight'		only neighbors within the mask
%|		'leak'		include neighbors just outside of mask (aspire)
%|		'none'		neither: user computes wk's elsewhere (default)
%|	'distance_power' val	1 classical (default), 2 possibly improved
%|	'order'	1 or 2		2 for 2nd order differences.  (default: 1)
%|
%| out
%|	ob	[np*noffset,np]		Fatrix object
%|
%| Copyright 2004-5-18, Jeff Fessler, University of Michigan

if nargin == 1 && streq(kappa, 'test'), Cdiff_test, return, end
if nargin < 1
	help(mfilename)
	if has_mex_jf, penalty_mex('help'), end
	error(mfilename)
end

%
% required input arguments
%
if numel(kappa) <= 3
	kappa = ones(kappa);
end

% option defaults
arg.offsets = [];
arg.order = 1;
arg.dosparse = 0;
arg.edge_type = 'none';
arg.distance_power = 1;
arg.Cpower = 1; % default of course is C^1

% parse optional name/value pairs
arg = vararg_pair(arg, varargin);

arg.dim_i = size(kappa);

%
% offsets to neighbors
%
if ischar(arg.offsets)
	if streq(arg.offsets, '2d,hvd')
		nx = arg.dim_i(1);
		arg.offsets = [1 nx nx+1 nx-1];
	else
		fail('unknown offset type %s', arg.offsets)
	end
end

if isempty(arg.offsets)
	if ndims(kappa) == 1
		arg.offsets = [1];
	elseif ndims(kappa) == 2
		nx = arg.dim_i(1);
		arg.offsets = [1 nx nx+1 nx-1]; % default 2D
	elseif ndims(kappa) == 3
		nx = arg.dim_i(1);
		ny = arg.dim_i(2);
		arg.offsets = [1 nx nx*ny]; % bare-bones 3D
	else
		error 'default for > 3D not done'
	end
end

arg.offsets = int32(arg.offsets);

%
% C_effective = Diag(scale) * C_raw
%
if streq(arg.edge_type, 'none')
	arg.scale = [];
else
	arg.scale = sprintf('wk,%s,%d', arg.edge_type, arg.order);
	arg.scale = penalty_mex(arg.scale, single(kappa), arg.offsets, ...
		arg.distance_power);
	arg.scale = sqrt(arg.scale);
end

arg.dim_o = [size(kappa) length(arg.offsets)];
arg.ndim = int32(ndims(kappa));

arg.mask = kappa ~= 0;
arg.np = sum(arg.mask(:));
arg.nd = prod(arg.dim_o);
ob = Fatrix([arg.nd arg.np], arg, 'caller', 'Cdiff', ...
	'forw', @Cdiff_forw, 'back', @Cdiff_back, 'power', @Cdiff_power);

if arg.dosparse
	ob = ob * speye(arg.np);
end


%
% Cdiff_power()
% for C.^2
function ob = Cdiff_power(ob, p)
ob.arg.Cpower = ob.arg.Cpower * p;

%
% Cdiff_forw()
% y = C * x
%
function y = Cdiff_forw(arg, x)

diff_str = sprintf('diff%d,forw%d', arg.order, arg.Cpower);

% if needed, convert concise column to array
flag_column = 0;
if size(x,1) == arg.np
	flag_column = 1;
	x = embed(x, arg.mask);
end

if ~isreal(x)
	yr = penalty_mex(diff_str, single(real(x)), arg.offsets, arg.ndim);
	yi = penalty_mex(diff_str, single(imag(x)), arg.offsets, arg.ndim);
	y = complex(yr, yi);
else
	y = penalty_mex(diff_str, single(x), arg.offsets, arg.ndim);
end
if flag_column % column(s) in yields column(s) out
	y = mult_rep(y, arg.scale);
	y = reshape(y, [arg.nd numel(y)/arg.nd]);
else
	if ~isempty(arg.scale)
		y = y .* arg.scale;
	end
end


%
% Cdiff_back()
% x = C' * y
%
function x = Cdiff_back(arg, y)

diff_str = sprintf('diff%d,back%d', arg.order, arg.Cpower);

flag_column = 0;
if size(y,1) == arg.nd
	flag_column = 1;
	y = reshaper(y, arg.dim_o);
end

y = mult_rep(y, arg.scale);
if ~isreal(y)
	xr = penalty_mex(diff_str, single(real(y)), arg.offsets, arg.ndim);
	xi = penalty_mex(diff_str, single(imag(y)), arg.offsets, arg.ndim);
	x = complex(xr, xi);
else
	x = penalty_mex(diff_str, single(y), arg.offsets, arg.ndim);
end

if flag_column
	x = reshape(x, [prod(arg.dim_i) numel(x)/prod(arg.dim_i)]);
	x = x(arg.mask,:);
end


%
% mult_rep()
% .* two arrays, expanding the second array if needed
%
function y = mult_rep(y, scale)
if isempty(scale)

elseif ndims(y) == ndims(scale)
	y = y .* scale;
else
	dim = size(y);
	dim([1:ndims(scale)]) = 1;
	scale = repmat(scale, dim);
	y = y .* scale;
end
