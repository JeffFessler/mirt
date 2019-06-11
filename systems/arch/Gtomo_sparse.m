 function ob = Gtomo_sparse(arg1, arg2, varargin)
%function ob = Gtomo_sparse(file(or matrix), mask, options)
%
% Construct Gtomo_sparse object, either from a .wtf file,
% or from a Matlab sparse matrix.
% This just uses an ordinary Matlab sparse matrix for the core!
% See Gtomo_sparse_test.m for example usage.
% You create an system object by calling:
%	G = Gtomo_sparse(file)
% and then you can use it thereafter by typing commands like
%	y = G * x.
% in
%	file	'string.wtf'		name of .wtf
%		or sparse matrix	in this case, must give the dimensions!
%	mask	[idim]			or empty to read from .wtf
%
% options
%	idim	[1,ndim_in]		input dimensions: nx,ny,nz etc.
%	odim	[1,ndim_out]		output dimensions: nb,na etc.
%	chat				verbosity
% out
%	ob	[nd,np]		nd = prod(odim), np = sum(mask(:))
%				so it is already "masked"
%
% If the input matrix is already "masked," then the mask argument is key.
%
% Copyright 2004-9-21, Jeff Fessler, The University of Michigan

if nargin == 1 & streq(arg1, ''), Gtomo_sparse_test, return, end

arg.mask = arg2;

% defaults
arg.chat = 0;
arg.idim = [];
arg.odim = [];
if nargin < 2, chat = 0; end

arg.blocks = {};	% place to store blocks of G for subset algorithms
arg = vararg_pair(arg, varargin);

%
% if input is a .wtf file
%
if ischar(arg1)
	arg.file = arg1;

	if ~isempty(arg.idim) | ~isempty(arg.odim)
		error 'idim / odim should not be given for .wtf'
	end
	[arg.G arg.idim(1) arg.idim(2) arg.odim(1) arg.odim(2)] = ...
		wtfmex('load', arg.file);

	% default mask from .wtf
	if isempty(arg.mask)
		tmp = sum(arg.G) > 0;
		arg.mask = reshape(tmp, arg.idim);
	end

	arg.G = arg.G(:,arg.mask(:));

%
% if input is a sparse matrix
%
elseif issparse(arg1)
	arg.G = arg1;
	if isempty(arg.idim) | isempty(arg.odim)
		error 'idim / odim must be given for sparse matrix!'
	end
	if prod(arg.odim) ~= size(arg.G,1)
		error 'bad row dimension'
	end

	% default mask is all
	if isempty(arg.mask)
		arg.mask = true(arg.idim);
	elseif ~islogical(arg.mask)
		error 'mask must be logical'
	end

	if length(arg.idim) ~= ndims(arg.mask) | any(arg.idim ~= size(arg.mask))
		disp(arg), error 'bad mask size'
	end

	if sum(arg.mask(:)) ~= size(arg.G,2)
		if size(arg.G,2) == numel(arg.mask)
			arg.G = arg.G(:, arg.mask(:)); % trick: compact size
		else
			disp(arg), error 'bad G size'
		end
	end

else
	error 'input must be filename or sparse matrix'
end

%
% build Fatrix object
%
arg.nd = prod(arg.odim);
arg.np = sum(arg.mask(:));
dim = [arg.nd arg.np];
ob = Fatrix(dim, arg, 'caller', mfilename, ...
        'forw', @Gtomo_sparse_forw, 'back', @Gtomo_sparse_back, ...
        'mtimes_block', @Gtomo_sparse_mtimes_block);


%
% Gtomo_sparse_forw(): y = G * x
%
function y = Gtomo_sparse_forw(arg, x)

% if needed, convert array to concise column
flag_array = 0;
if size(x,1) ~= arg.np
	flag_array = 1;
	x = reshape(x, numel(arg.mask), []);
	x = x(arg.mask(:),:); % [np, L]
end

y = arg.G * x;

if flag_array
	y = reshaper(y, arg.odim);
end


%
% Gtomo_sparse_back(): x = G' * y
%
function x = Gtomo_sparse_back(arg, y)

flag_array = 0;
if size(y,1) ~= arg.nd
	flag_array = 1;
	y = reshape(y, arg.nd, []);
end

x = arg.G' * y;

if flag_array
	x = embed(x, arg.mask);
end


%
% Gtomo_sparse_mtimes_block()
%
function y = Gtomo_sparse_mtimes_block(arg, is_transpose, x, istart, nblock)
error 'not done'
