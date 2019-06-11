 function ob = Cdiffs(isize, varargin)
%function C1 = Cdiffs(isize, [options])
%|
%| Construct C1 object that can compute C1 * x and the adjoint C1' * d
%| for a "finite differences" matrix C for roughness penalty regularization.
%| This "stacks up" multiple Cdiff1() objects, e.g., akin to using vertcat(),
%| for possible internal use by roughness penalty objects like Reg1().
%|
%| One can use this for "bilateral total variation (TV)" regularization.
%|
%| Caution: for large problem sizes, computing C1' * (C1 * x) will require
%| #offsets * #pixels intermediate memory to store C1*x, which may be too much.
%| To save memory, use Reg1 which can compute \sum_{m=1}^M C1_m' * (C1_m * x).
%| This Cdiffs object is provided mostly for completeness and for small cases.
%|
%| Caution: for most types the "differences" here wrap around image edges.
%| User must combine with Rweights to multiply such differences by zero.
%| Use C = Gdiag(wt) * Cdiffs(...).  See inside Reg1 for examples.
%|
%| in
%|	isize	[]	vector of object dimensions (N), e.g., [64 64]
%|			if empty, infer from 'mask' option
%|
%| options
%|	'type_diff'	see Cdiff1.m (default: '' defers to Cdiff1)
%|	'offsets' [M]	offsets to "M" neighbors; see penalty_offsets()
%|	'order'	1 or 2	1st- or 2nd-order differences.  (default: 1)
%|	'mask'	[(N)]	logical support
%|	'class'	''	'fatrix2' or 'Fatrix' or '' (default: defer to Cdiff1)
%|	'odim_squeeze' 1|0	if M=1 then use odim = [(N)] not [(N) M]
%|				default: 1
%|
%| out
%|	C1	[*N * M, np]	fatrix2 or Fatrix object; np = sum(mask(:))
%|				also works on arrays: [(N) (L)] -> [(N) M (L)]
%|				or sparse matrix for type_diff == 'spmat'
%|			(sparse version handles only vectors, not arrays)
%|
%| Copyright 2006-12-4, Jeff Fessler, University of Michigan

if nargin == 1 && streq(isize, 'test'), Cdiffs_test, return, end
if nargin < 1, ir_usage, end
%if has_mex_jf, penalty_mex('help'), end

% option defaults
arg.type_diff = '';
arg.offsets = [];
arg.mask = [];
arg.order = 1;
arg.odim_squeeze = true;
arg.class = ''; % defer to Cdiff1

% parse optional name/value pairs
arg = vararg_pair(arg, varargin);

if isempty(isize)
	if isempty(arg.mask)
		fail 'must provide at least one of isize or ''mask'''
	else
		isize = size(arg.mask);
		if numel(isize) == 2 && isize(end) == 1
			warn 'inferring isize from 1D mask is ambiguous'
		end
	end
elseif ~isempty(arg.mask)
	msize = size(arg.mask);
	if numel(isize) == 1 && numel(msize) == 2 ...
		&& msize(1) == isize && msize(2) == 1 % 1D special case
		; % ok
	elseif ~isequal(msize, isize)
		fail 'size(mask) must match the ''isize'' argument'
	end
end

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
	if ~isempty(arg.mask)
		ob = ob(:,arg.mask(:));
	end

% typical object case
else
	arg.isize = isize;
	arg.Cc = cell(MM,1);
	for mm=1:MM
		arg.Cc{mm} = Cdiff1(isize, 'type_diff', arg.type_diff, ...
			'class', arg.class, ...
			'offset', arg.offsets(mm), 'order', arg.order);
	end

	if isempty(arg.mask)
		arg.mask = true([isize 1]);
	end

	if MM >= 1
		arg.class = class(arg.Cc{1}); % trick: defer to Cdiff1
	else % for null case
%		arg.class = 'fatrix2'; % todo: this fails for Reg1 test
		arg.class = 'Fatrix';
	end

	switch arg.class
	case 'fatrix2' % cannot use vertcat because Cdiff1 lacks mask
		if arg.odim_squeeze && MM == 1
			odim = [arg.isize];
		else
			odim = [arg.isize MM];
		end
		ob = fatrix2('arg', arg, 'odim', odim, ...
			'forw', @Cdiffs_forw, ...
			'back', @Cdiffs_back, ...
			'power', @Cdiffs_power, 'abs', @Cdiffs_abs);
	case 'Fatrix'
		arg.np = sum(arg.mask(:));

		dim = [prod(isize)*MM arg.np];
		ob = Fatrix(dim, arg, 'caller', 'Cdiffs', ...
			'forw', @Cdiffs_forw_Fatrix, ...
			'back', @Cdiffs_back_Fatrix, ...
			'power', @Cdiffs_power, 'abs', @Cdiffs_abs);
	otherwise
		fail('unknown class "%s"', arg.class)
	end

%	ob = block_fatrix(ob, 'type', 'col'); % old approach
end


%
% Cdiffs_forw(): y = A * x
%
function y = Cdiffs_forw(arg, x)

MM = numel(arg.Cc);
yy = cell(MM,1);
for mm=1:MM
	yy{mm} = arg.Cc{mm} * x;
end
dim_cat = numel(arg.isize) + 1;
y = cat(dim_cat, yy{:});


%
% Cdiffs_forw_Fatrix(): y = A * x
%
function y = Cdiffs_forw_Fatrix(arg, x)

[x ei] = embed_in(x, arg.mask, arg.np); % [(N) *L]
LL = size(x, 1+length(arg.isize));

MM = length(arg.Cc);

y = zeros([prod(arg.isize)*LL MM]); % [*N * *L, M]
for mm=1:MM
	tmp = arg.Cc{mm} * x; % [(N) *L]
	y(:,mm) = tmp(:);
end

if LL > 1
	y = reshape(y, [prod(arg.isize) LL MM]); % [*N *L M]
	y = permute(y, [1 3 2]); % [*N M *L]
end
y = reshape(y, [arg.isize MM LL]); % [(N) M *L]

y = ei.shape(y); % [*N * M, (L)] or [(N) M (L)]


%
% Cdiffs_back(): x = A' * y
%
function x = Cdiffs_back(arg, y)

MM = length(arg.Cc);
y = reshapee(y, [], MM); % [*N M]

x = 0;
for mm=1:MM
	tmp = reshape(y(:,mm), [arg.isize 1]); % [(N)]
	tmp = arg.Cc{mm}' * tmp; % [(N)]
	tmp = tmp .* arg.mask;
	x = x + tmp;
end


%
% Cdiffs_back_Fatrix(): x = A' * y
%
function x = Cdiffs_back_Fatrix(arg, y)

MM = length(arg.Cc);
[y eo] = embed_out(y, [arg.isize, MM]); % [(N) M *L]
LL = size(y, 2+length(arg.isize));

y = reshape(y, [prod(arg.isize) MM LL]); % [*N M *L]

if LL > 1
	y = permute(y, [1 3 2]); % [*N *L M]
	y = reshape(y, [prod(arg.isize)*LL MM]); % [*N * *L, M]
end

x = 0;
for mm=1:MM
	tmp = reshape(y(:,mm), [arg.isize LL]); % [(N) *L]
	tmp = arg.Cc{mm}' * tmp; % [(N) *L]
	x = x + tmp;
end

if LL > 1
	x = x .* repmat(arg.mask, [ones(1,ndims(arg.mask)) LL]);
else
	x = x .* arg.mask;
end
x = eo.shape(x, arg.mask, arg.np); % [*N (L)] or [(N) (L)]


%
% Cdiffs_abs()
%
function ob = Cdiffs_abs(ob)
MM = length(ob.arg.Cc);
Ca = cell(MM,1);
for mm=1:MM
	Ca{mm} = abs(ob.arg.Cc{mm});
end
ob.arg.Cc = Ca;


%
% Cdiffs_power()
% for C.^2
%
function ob = Cdiffs_power(ob, p)
MM = length(ob.arg.Cc);
Cp = cell(MM,1);
for mm=1:MM
	Cp{mm} = ob.arg.Cc{mm} .^ p;
end
ob.arg.Cc = Cp;


%
% Cdiffs_test()
%
function Cdiffs_test
ig = image_geom('nx', 8, 'ny', 6, 'dx', 1);
%ig.mask = ig.circ > 0;
ig.mask(3) = 0; % stress
%im(ig.mask)

% x = ig.unitv;
rng(0)
x = rand(ig.dim);

list_class = {'fatrix2', 'Fatrix'};

for ic = 1:numel(list_class)
 for order=1:2
	pr order
	args = {ig.dim, 'order', order, 'mask', ig.mask, ...
		'class', list_class{ic}, 'type_diff'};
	Ci = Cdiffs(args{:}, 'ind'); % a basic one

	wt = Rweights(ig.mask, Ci.arg.offsets, 'type_wt', 'array', ...
		'order', order, 'distance_power', 0); % needed for comparisons
	wt = ig.shape(wt); 

	y1 = Ci * x;
	y1w = wt .* y1;
	z1 = Ci' * y1;
	z1w = Ci' * y1w;

%{
	switch list_class{ic}
	case 'Fatrix'
		Fatrix_test_basic(Ci, ig.mask)
	case 'fatrix2'
		fatrix2_tests(Ci)
	end
%}

	Ci = Ci(:,:);

	types = Cdiff1('types'); % all of them except spmat
%		{'def', 'ind', 'mex', 'sparse', 'convn', 'imfilter'};

	for it=1:numel(types)
%		if streq(types, 'spmat'), continue, end
		C = Cdiffs(args{:}, types{it});
		y2 = C * x;
		z2 = C' * y1;
		z2w = C' * y1w;

		switch types{it}
		case {'circshift', 'convn', 'imfilter'} % trick: require wt
			equivs(y1w, wt .* y2)
			equivs(z1w, z2w)
			jf_equal(diag(wt(:)) * C(:,:), diag(wt(:)) * Ci)
		otherwise
			equivs(y1, y2)
			equivs(z1, z2)
			Ct = C(:,:);
			jf_equal(Ct, Ci)
		end

		switch list_class{ic}
		case 'Fatrix'
			Fatrix_test_basic(C, ig.mask)
		case 'fatrix2'
			fatrix2_tests(C)
		end

		Cdiff1_test1(C) % abs, squared, adjoint
	end

	if 0
		im plc 1 2
		im(1, y1)
		im(2, z1)
	end

	% test sparse matrix too
	Cs = Cdiffs(args{:}, 'sparse');
	Cz = Cdiffs(args{:}, 'spmat');
	Cf = Cs(:,:);
	jf_equal(Cf, Cz)

	% abs
	Ca = abs(Cs);
	Cf = Ca(:,:);
	jf_equal(Cf, abs(Cz))

	% abs
	Cp = Cs .^ 2;
	Cf = Cp(:,:);
	jf_equal(Cf, Cz .^ 2)
 end
end
