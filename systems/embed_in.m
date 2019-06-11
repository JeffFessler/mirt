  function [x, ei_fun] = embed_in(x, mask, np)
%|function [x, ei_fun] = embed_in(x, mask, np)
%|
%| The matrix-vector versions of the linear forward models used herein
%| expect an input x of size [np *L] and produce an output y of size [nd *L],
%| where np = sum(mask(:)) is the number of estimated pixel values,
%| and nd is the number of data points.
%| However, for convenience it is also desirable for the overloaded
%| mtimes operation to map a [(N) (L)] input into a [(M) (L)] output,
%| where (N) is the d_in -dimensional input size
%| and (M) is the d_out -dimensional output size,
%| and (L) reflects the possibility of multiple inputs, e.g.: A * [u v w].
%|
%| Furthermore, most of the code for forward models is designed to
%| map input [(N) *L] to output [(M) *L].
%|
%| in
%|	x	[np (L)] or [(N) (L)]	input image(s), possibly as columns
%|	mask	[(N)]			logical support array
%|	np				sum(mask(:))
%| out
%|	x	[(N) *L]		output images, as array
%|	ei_fun	strum object with methods:
%|		y = ei_fun.shape(y)	reshape y from [(M) *L] to be either
%|					[nd (L)] or [(M) (L)], depending on x
%|
%| Copyright 2006-12-9, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), embed_in_test, clear x, return, end
if nargin < 2, help(mfilename), error(mfilename), end
if nargin < 3, np = sum(mask(:)); end

% convert input to [(N) *L]
state.column = false;
if size(x,1) == np % convert [np (L)] to [(N) *L]
	state.column = true;
	dimx = size(x);
	state.diml = dimx(2:end); % [(L)]
	x = embed(x, mask, '*dim', 0); % [(N) *L]

else % convert [(N) (L)] to [(N) *L]
	dimi = size(x);
	if length(dimi) < ndims(mask), fail('bad image size'), end

	size1_mask = size(mask);
	ndims1_mask = ndims(mask);
	if size1_mask(end) == 1 % trick: handle '1d' mask well
		ndims1_mask = ndims1_mask - 1;
		size1_mask = size1_mask(1:end-1);
	end
	jf_equal(dimi(1:ndims1_mask), size1_mask)
%	jf_equal(dimi(1:length(size(mask))), size(mask)) % pre 2012-06-03
	state.diml = dimi((ndims1_mask+1):end); % (L)
%	x = reshape(x, [size(mask) prod(state.diml) 1]); % [(N) *L] % pre
	x = reshape(x, [size1_mask prod(state.diml) 1]); % [(N) *L]

%	if ~dims_same(x, mask, 'up_to_dim', ndims(mask))
%		error(['dimension mismatch.  x: ' num2str(size(x), ' %0d') ...
%			', mask: ' num2str(size(mask), ' %0d')])
%	end
end

ei_fun = strum(state, {'shape', @embed_in_shape, '(y)'});


%
% embed_in_shape()
%
function y = embed_in_shape(state, y)

diml = state.diml;

if state.column % column in yields column out, i.e., [(M) *L] to [*M (L)]
	if any(diml > 1)
		diml = num2cell(diml);
		y = reshape(y, [], diml{:}); % [*M (L)]
	else
		y = y(:); % [*M,1]
	end

else % [(M) *L] to [(M) (L)]
	if any(diml > 1)
		dimy = size(y);
		if dimy(end) ~= prod(diml), error 'diml bug', end
		dimm = dimy(1:end-1); % (M)
		y = reshape(y, [dimm diml]);
	end
end


function embed_in_test
ig = image_geom('nx', 10, 'ny', 8, 'dx', 1);
ig.mask = ig.circ > 0;

x = ig.unitv;
dl = [2 3];
c = repmat(x(ig.mask), [1 dl]); % [np (L)]
[x1 ei] = embed_in(c, ig.mask);
t = ei.shape(ones(4, 5, prod(dl)));
jf_equal(size(t), [4*5 dl])

x2 = repmat(x, [1 1 dl]); % [(N) (L)]
%size(x2)
[x3 ei] = embed_in(x2, ig.mask);
jf_equal(size(x3), [ig.dim prod(dl)])
t = ei.shape(ones(4, 5, prod(dl)));
jf_equal(size(t), [4 5 dl])
