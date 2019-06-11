 function y = embed_mult(fun, arg, is_transpose, x, istart, nblock, ...
	mask, np, odim, is_array, varargin)
%function y = embed_mult(fun, arg, is_transpose, x, istart, nblock, ...
%|	mask, np, odim, is_array, varargin)
%|
%| The matrix-vector versions of the linear forward models used herein
%| expect an input x of size [np *L] and produce an output y of size [nd *L],
%| where np = sum(mask(:)) is the number of estimated pixel values.
%| However, for convenience it is also desirable for the overloaded
%| mtimes operation to map a [(N) (L)] input into a [(M) (L)] output,
%| where (N) is the d_in - dimensional input size
%| and (M) is the d_out - dimensional output size,
%| and (L) reflects the possibility of multiple inputs, e.g.: A * [u v w].
%|
%| Furthermore, most of the code for forward models is designed to
%| map input [(N) *L] to output [(M) *L].
%|
%| This is a generic interface to work with these possibilities.
%|
%| in
%|	fun		function handle, call:
%|			y = fun(arg, is_tranpose, x, istart, nblock, varargin)
%|	arg	struct
%|	x	[]	input image
%|	is_transpose	0|1
%|	istart,nblock	for block operations
%|	mask	[(N)]	logical support array
%|	np		sum(mask(:))
%|	odim		(L)
%|	is_array	1 for now
%| out
%|	y	[]	output data
%|
%| Copyright 2006-12-9, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fun, 'test'), embed_mult_test, return, end
if nargin < 10, help(mfilename), error(mfilename), end

switch is_array
case 1
	switch is_transpose
	case 0
		y = embed_mult_forw_array(fun, arg, x, mask, np, ...
			istart, nblock, varargin{:});
	case 1
		y = embed_mult_back_array(fun, arg, x, mask, odim, ...
			istart, nblock, varargin{:});
	otherwise
		error 'bug'
	end

otherwise
	error 'only array done'
end


%
% embed_mult_forw_array()
% For the case that the function maps input [(N),*L] to output [(M),*L].
%
function y = embed_mult_forw_array(fun, arg, x, mask, np, istart, nblock, varargin)

% convert input to [(N) *L]
flag_column = 0;
if size(x,1) == np % convert [np (L)] to [(N) *L]
	flag_column = 1;
	dimx = size(x);
	diml = dimx(2:end); % [(L)]
	x = embed(x, mask, '*dim', 0); % [(N) *L]

else % convert [(N) (L)] to [(N) *L]
	dimi = size(x);
	diml = dimi(ndims(mask)+1:end); % (L)
	x = reshape(x, [size(mask) prod(diml) 1]); % [(N) *L]

%	if ~dims_same(x, mask, 'up_to_dim', ndims(mask))
%		error(['dimension mismatch.  x: ' num2str(size(x), ' %0d') ...
%			', mask: ' num2str(size(mask), ' %0d')])
%	end
end

y = fun(arg, 0, x, istart, nblock, varargin{:}); % [(M) *L]

%na = arg.nn(end); % last index must be the one over which we subsetize
%ia = istart:nblock:na;
%y = reshape(y, prod(arg.nn(1:end-1)), na, []);
%y = y(:,ia,:);

if flag_column % column in yields column out, i.e., [(M) *L] to [*M (L)]
	if any(diml > 1)
		diml = num2cell(diml);
		y = reshape(y, [], diml{:}); % [*M (L)]
	else
		y = y(:); % [*M 1]
	end
else % [(M) *L] to [(M) (L)]
	if any(diml > 1)
		dimy = size(y);
		if dimy(end) ~= prod(diml), error 'diml bug', end
		dimm = dimy(1:end-1); % (M)
		y = reshape(y, [dimm diml]);
	end
end


%
% embed_mult_back_array()
% Adjoint when the forward function maps input [(N) *L] to output [(M) *L].
%
function x = embed_mult_back_array(fun, arg, y, mask, odim, istart, nblock, varargin)

% convert input to [(M) *L]
if nblock == 1
	flag_column = 0;
	if size(y,1) == prod(odim) % [*M (L)]
		flag_column = 1;
		diml = size(y); diml = diml(2:end); % (L)
	else % [(M) (L)]
		diml = size(y); diml = diml(length(odim)+1:end); % (L)
	end
	y = reshape(y, [odim prod(diml)]); % [(M) *L]

else % y is [(Mi,Ms) (L)]
	na = odim(end); % Ms: last index must be the one over which we subsetize
	ia = istart:nblock:na;
	n1 = prod(odim(1:end-1)); % *Mi

	flag_column = 0;
	if size(y,1) == n1 * length(ia) % [*Mi * Ms, (L)]
		flag_column = 1;
		diml = size(y); diml = diml(2:end); % (L)
	else % [(Mi,Ms),(L)]
		diml = size(y); diml = diml(length(odim)+1:end); % (L)
	end

	y = reshape(y, n1, length(ia), []); % [*Mi Ms *L]
	yy = zeros([n1 na prod(diml)], class(y)); % [*Mi Ms *L]
	yy(:,ia,:) = y;
	dim = num2cell(odim);
	y = reshape(yy, dim{:}, []); % [(M) *L]
end

x = fun(arg, 1, y, istart, nblock, varargin{:}); % [(N) *L]

if flag_column
	x = reshape(x, numel(mask), []); % [*N *L]
	x = x(mask,:); % [np *L]
	x = reshape(x, [size(x,1) diml]); % [np (L)]
else
	x = reshape(x, [size(mask) diml]); % [(N) (L)]
end


function embed_mult_test
Gtomo2_dscmex test
