 function ob = Cdiffs(isize, varargin)
%function C1 = Cdiffs(isize, [options])
%
% Construct C1 object that can compute C1 * x and the adjoint C1' * d
% for a "differencing" matrix C for roughness penalty regularization.
% This object "stacks up" multiple Cdiff1() objects using block_fatrix()
% for possible internal use by the roughness penalty objects.
%
% Caution: for large problem sizes, computing C1' * (C1 * x) will require
% #offsets * #pixels intermediate memory to store C1*x, which may be too much.
% Instead, one should compute \sum_{m=1}^M C1_m' * (C1_m * x).
% So this is really mostly for completeness.
%
% in
%	isize	[]		vector of object dimensions (N), e.g., [64 64]
%
% options
%	'type_diff'		see Cdiff1.m
%	'offsets' [M]		offsets to "M" neighbors; see penalty_offsets()
%	'order'	1 or 2		1st- or 2nd-order differences.  (default: 1)
%
% out
%	C1	[*N * M;*N]	Fatrix object.  trick: also works [(N)*M;(N)] ??
%				or sparse matrix for type_diff == 'spmat'
%
% Copyright 2006-12-4, Jeff Fessler, The University of Michigan

if nargin == 1 & streq(isize, 'test'), Cdiffs_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end
%if has_mex_jf, penalty_mex('help'), end

% option defaults
arg.type_diff = '';
arg.offsets = [];
arg.order = 1;

% parse optional name/value pairs
arg = vararg_pair(arg, varargin);

% offsets to neighbors
arg.offsets = penalty_offsets(arg.offsets, isize);

MM = length(arg.offsets);

% sparse matrix case
if streq(arg.type_diff, 'spmat')
	ob = [];
	for mm=1:MM
		ob = [ob; Cdiff1(isize, 'type_diff', arg.type_diff, ...
			'offset', arg.offsets(mm), 'order', arg.order)];
	end

% typical object case
else
	ob = cell(MM,1);
	for mm=1:MM
		ob{mm} = Cdiff1(isize, 'type_diff', arg.type_diff, ...
			'offset', arg.offsets(mm), 'order', arg.order);
	end

	ob = block_fatrix(ob, 'type', 'col');
end


function Cdiffs_test
ig = image_geom('nx', 8, 'ny', 6, 'dx', 1);
types = {'def', 'ind', 'mex', 'sparse'};
for it=1:length(types)
	C = Cdiffs(ig.dim, 'type_diff', types{it}, 'order', 1);
end

% test sparse matrix too
Cz = Cdiffs(ig.dim, 'type_diff', 'spmat', 'order', 1);
C = C(:,:);
if ~isequal(C, Cz), error 'bug', end

im plc 1 2
y = C * ig.unitv(:);
im(1, ig.shape(y))
z = C' * y;
z = ig.embed(z);
im(2, z)
