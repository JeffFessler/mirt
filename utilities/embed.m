 function ff = embed(x, mask, varargin)
%function ff = embed(x, mask, varargin)
%| embed x in nonzero elements of (logical) mask
%| in
%|	x	[np (L)]	the "nonzero" pixels (lexicographically stacked)
%|	mask	[(Nd)]		logical array, np = sum(mask)
%| option
%|	'*dim'	{0|1}		0: [(N) (L)] (default); 1: return [(N) *L]
%| out
%|	ff	[(Nd) (L)]	viewable image(s)
%|				if input is sparse, output is full double
%|
%| See also: masker
%|
%| Copyright 2000-9-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), embed_test, return, end
if nargin < 2, ir_usage, end

if isempty(mask) % trick: handle empty mask quickly, ignoring '*dim'
	ff = x;
return
end

arg.prod_dim = 0;
if nargin > 2
	arg = vararg_pair(arg, varargin, 'subs', {'*dim', 'prod_dim'});
end

if ~islogical(mask), error 'mask must be logical', end

dimx = size(x);
cl = class(x);
if issparse(x)
	cl = 'double';
end

pL = prod(dimx(2:end));

if islogical(x)
%	cl = 'double';
	ff = false([numel(mask) pL]); % [np *L]
else
	ff = zeros([numel(mask) pL], cl); % [np *L]
end

%if is_pre_v7
%	ff = zeros([numel(mask) pL]); % [np *L]
%else
%	ff = zeros([numel(mask) pL], cl); % [np *L]
%end

if pL > 1
	ff(mask(:),:) = reshapee(x, [], pL);
else
	ff(mask(:),1) = x;
end

if ~arg.prod_dim
	if ndims(mask) == 2 && size(mask,2) == 1 % 1d cases, 2008-12-14
		% trick: omit '1' in 2nd dimension for 1d cases:
		ff = reshape(ff, [size(mask,1) dimx(2:end)]); % [(Nd) (L)]
	else
		ff = reshape(ff, [size(mask) dimx(2:end)]); % [(Nd) (L)]
	end
end


function embed_test
ig = image_geom('nx', 512, 'ny', 500, 'dx', 1);
ig.mask = ig.circ > 0;
ig.mask = conv2(double(ig.mask), ones(2), 'same') > 0;

x = [1:sum(ig.mask(:))]';
cpu etic
f1 = embed(x, ig.mask);
cpu etoc 'time for one'
f2 = embed([x x], ig.mask);
cpu etoc 'time for two'
jf_equal(f1, f2(:,:,2))

x = repmat(x, [1 2 3]);
f3 = embed(x, ig.mask);
jf_equal(f1, f3(:,:,2))

x = ig.unitv;
x = x(ig.mask);
f4 = ig.embed(x);
x = sparse(double(x));
% f5 = ig.embed(x); % no! trick in image_geom.m for sparse x
f5 = embed(x, ig.mask);
jf_equal(f4, f5)

% test 1d case
mask = logical([0 1 1 1 0]');
x = [1 2 3]';
jf_equal(embed(x,mask), [0 1 2 3 0]')
x = [1 2 3]' * [1 1]; % [np 2]
jf_equal(embed(x,mask), [0 1 2 3 0]'*[1 1])
