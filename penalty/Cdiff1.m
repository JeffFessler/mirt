 function ob = Cdiff1(isize, varargin)
%function C1 = Cdiff1(isize, [options])
%|
%| Construct Cdiff1 object that can compute C1 * x and the adjoint C1' * d
%| for a "finite differencing" matrix C1 for roughness penalty regularization.
%| This object provides the finite difference in only a *single* direction,
%| so it is designed to be used internally by roughness penalty objects.
%| 
%| Cdiffs "stacks these up" which is more useful (if memory permits).
%|
%| For the usual 1st-order differences, each row of C1 is either all zeros
%| (near border) or all zeros except for a single +1 and -1 value.
%|
%| Caution: most of these disregard image borders so users will probably
%| want to combine this with Rweights().  This is done already in Reg1().
%| Also, this does not support any 'mask' option.  Use Cdiffs() for that.
%|
%| in
%|	isize	[]		vector of object dimensions (N), e.g., [64 64]
%|
%| options
%|	'type_diff'		method for performing differencing:
%|				'def' | '' - default : try in this order:
%|					'mex' 'circshift'
%|				'circshift' - using matlab's circshift()
%|				'diff' - using matlab diff() todo - not done!
%|				'ind' - using matlab indexing
%|				'for1' - using matlab "for" loop (slow)
%|				'mex' - using penalty_mex file (preferable)
%|				'sparse' - using internal Gsparse object
%|				'spmat' - return sparse matrix itself
%|				'imfilter' - use matlab's imfilter()
%|					(requires image processing toolbox)
%|				'convn' - use matlab's convn, e.g. with [-1 1]
%|					(very slow - do not use!)
%|	'offset' [int]		neighbor offset. (default: [1], left neighbor)
%|				Scalar, or displacement [dx dy dz ...] with length(isize).
%|				If 0, then C1 = Identity.
%|	'order'	1 or 2		1st- or 2nd-order differences.  (default: 1)
%|	'class'			'fatrix2' (default) or 'Fatrix'
%|
%| out
%|	C1	[*N *N]		fatrix2 or Fatrix object
%|				or sparse matrix, if 'spmat' option.
%|				trick: also works [(N) (L)] -> [(N) (L)]
%|					or [*N (L)] -> [*N (L)]
%|
%| The name "C1" reminds us the elements are +/- 1 and 0 (for order=1).
%|
%| For large images (where it matters) the compute times are roughly
%|	convn > for1 > ind > imfilter > sparse > circshift > mex
%| So it seems that only the 'mex' and 'circshift' options are useful.
%| However, it may be cpu dependent so use "Cdiff1_tune" to compare.
%|
%| Copyright 2006-11-29, Jeff Fessler, University of Michigan

if nargin == 1 && streq(isize, 'test'), Cdiff1_test, return, end
if nargin < 1, ir_usage, end

if streq(isize, 'types') % list from fastest to slowest order (on 2008 mac pro)
	ob = {'mex', 'circshift', 'imfilter', 'sparse', 'for1', 'ind', 'convn'};
	return
end

% required input argument(s)
arg.isize = isize;

% option defaults
arg.type_diff = 'def';
arg.offset = [1];
arg.order = 1;
arg.class = ''; % see below

% parse options
arg = vararg_pair(arg, varargin);

if isempty(arg.class)
	arg.class = 'fatrix2';
end

if length(arg.offset) > 1 % displacement vector instead of offset scalar
	if ~isequal(size(arg.offset), size(arg.isize)), error 'offset size', end
	arg.offset = sum(cumprod([1 arg.isize(1:end-1)]) .* arg.offset);
end
arg.displace = penalty_displace(arg.offset, arg.isize);

arg.is_abs = false;
arg.Cpower = 1; % start with C = C^1
arg.nn = prod(arg.isize);
arg.dim = [arg.nn arg.nn];

if isempty(arg.type_diff) || streq(arg.type_diff, 'def') % default type
	if exist('penalty_mex') == 3
		arg.type_diff = 'mex';
%	elseif ir_has_imfilter
%		arg.type_diff = 'imfilter';
	else
		arg.type_diff = 'circshift';
	end
end

% make sure mex can be supported if requested.  if not, revert gracefully
if streq(arg.type_diff, 'mex') && exist('penalty_mex') ~= 3
	warn 'no penalty_mex so ''mex'' unavailable; using ''circshift'''
	arg.type_diff = 'circshift';
end

if streq(arg.type_diff, 'imfilter') && ~ir_has_imfilter
	warn 'no imfilter (no image processing toolbox); using ''circshift'''
	arg.type_diff = 'circshift';
end

% identity matrix is special case
if arg.offset == 0
	if arg.order ~= 0, warn 'only order=0 for identity', end

	switch arg.type_diff
	case 'spmat'
%		n = arg.nn;
%		ob = sparse(1:n, 1:n, ones(n,1), n, n);
		ob = speye(arg.nn);
	case 'sparse'
		ob = Gdiag(ones([arg.isize 1]), 'class', arg.class);
%		ob = diag_sp(ones(arg.nn,1));
	otherwise
		switch arg.class
		case 'fatrix2' % identity
			handle_forw = @(arg,x) x;
			handle_back = @(arg,y) y;
			handle_abs = @(ob) ob;
			handle_pow = @(ob, sup) ob;
			ob = fatrix2('idim', arg.isize, 'odim', arg.isize, ...
				'does_many', true, 'caller', 'Cdiff1:ident', ...
				'forw', handle_forw, 'back', handle_back, ...
				'abs', handle_abs, 'power', handle_pow);
		case 'Fatrix'
			ob = Fatrix(arg.dim, arg, 'caller', 'Cdiff1:ident', ...
				'forw', @Cdiff1_ident_dup, ...
				'back', @Cdiff1_ident_dup, ...
				'abs', @Cdiff1_abs, 'power', @Cdiff1_power);
		otherwise
			fail('class %s unknown', arg.class)
		end
	end
return
end

if ~any(arg.order == [1 2])
	error 'only 1st and 2nd order done (for nonzero offsets)'
end

caller = ['Cdiff1:' arg.type_diff];
Cdiff1_Fatrix = @(arg, forw, back) ...
	Fatrix(arg.dim, arg, 'caller', caller, ...
		'forw', forw, 'back', back, ...
		'power', @Cdiff1_power, 'abs', @Cdiff1_abs);
Cdiff1_fatrix2 = @(arg, does_many, forw, back) ...
	fatrix2('arg', arg, 'caller', caller, 'does_many', does_many, ...
		'idim', arg.isize, 'odim', arg.isize, ...
		'forw', forw, 'back', back, ...
		'power', @Cdiff1_power, 'abs', @Cdiff1_abs);

switch arg.type_diff

case 'circshift'
	arg = Cdiff1_coef_setup(arg);
	switch arg.class
	case 'fatrix2'
		ob = Cdiff1_fatrix2(arg, true, ...
			@Cdiff1_cs_forw, @Cdiff1_cs_back);
	case 'Fatrix'
		ob = Cdiff1_Fatrix(arg, ...
			@Cdiff1_cs_forw_Fatrix, @Cdiff1_cs_back_Fatrix);
	otherwise, fail bug
	end

case {'convn', 'imfilter'}
	arg = Cdiff1_coef_setup(arg);
	arg = Cdiff1_filt_setup(arg);
	switch arg.class
	case 'fatrix2'
		forw = @(arg, x) arg.filt_forw(x, ...
			Cdiff1_coef_powabs(arg.coef, arg));
		back = @(arg, y) arg.filt_forw(y, ... % trick!
			flipdims(Cdiff1_coef_powabs(arg.coef, arg)));
		ob = Cdiff1_fatrix2(arg, true, forw, back);
	case 'Fatrix'
		ob = Cdiff1_Fatrix(arg, ...
			@Cdiff1_filt_forw_Fatrix, @Cdiff1_filt_back_Fatrix);
	otherwise, fail bug
	end

case 'diff'
	if arg.order ~= 1 || ~streq(arg.class, 'fatrix2')
		fail 'only order=1 and fatrix2 for "diff"'
		% codo: could support 2nd order using diff...
	end
	displace = arg.displace; % should be all zeros except for a single "1"
	arg.diff_dim = find(displace == 1);
	displace(arg.diff_dim) = 0;
	if numel(arg.diff_dim) ~= 1 || any(displace)
		pr arg.displace
		fail('"diff" supports only simple offsets along coordinates')
	end

	ob = Cdiff1_fatrix2(arg, true, @Cdiff1_diff_forw, @Cdiff1_diff_back);

case 'for1' % trick: just using "ind" for back because for1 is slow anyway
	switch arg.class
	case 'fatrix2'
		forw = @(arg, x) Cdiff1_for1_forw(arg, x);
		back = @(arg, y) Cdiff1_ind_back(arg, y); % !
		ob = Cdiff1_fatrix2(arg, false, forw, back);
	case 'Fatrix'
		ob = Cdiff1_Fatrix(arg, ...
			@Cdiff1_for1_forw_Fatrix, @Cdiff1_ind_back_Fatrix); % !
	otherwise, fail bug
	end

case 'ind'
	switch arg.class
	case 'fatrix2'
		forw = @(arg, x) Cdiff1_ind_forw(arg, x);
		back = @(arg, y) Cdiff1_ind_back(arg, y);
		ob = Cdiff1_fatrix2(arg, false, forw, back);
	case 'Fatrix'
		ob = Cdiff1_Fatrix(arg, ...
			@Cdiff1_ind_forw_Fatrix, @Cdiff1_ind_back_Fatrix);
	otherwise, fail bug
	end

case 'mex'
	arg.isize32 = int32(length(arg.isize));
	arg.offset = int32(arg.offset);
	switch arg.class
	case 'fatrix2'
		forw = @(arg, x) Cdiff1_mex_forw(arg, x);
		back = @(arg, y) Cdiff1_mex_back(arg, y);
		ob = Cdiff1_fatrix2(arg, true, forw, back);
	case 'Fatrix'
		ob = Cdiff1_Fatrix(arg, ...
			@Cdiff1_mex_forw_Fatrix, @Cdiff1_mex_back_Fatrix);
	otherwise, fail bug
	end

case 'sparse'
	arg.C = Cdiff1_sp_setup(arg);
	switch arg.class
	case 'fatrix2'
		ob = Gmatrix(arg.C, 'idim', arg.isize, 'odim', arg.isize);
	case 'Fatrix'
		arg.C = Gsparse(arg.C, 'idim', arg.isize, 'odim', arg.isize);
		ob = Cdiff1_Fatrix(arg, @Cdiff1_sp_forw, @Cdiff1_sp_back);
	otherwise, fail bug
	end

case 'spmat'
	ob = Cdiff1_sp_setup(arg);

otherwise
	fail('bad type %s', arg.type_diff)
end


% Cdiff1_abs()
% for abs(C)
function ob = Cdiff1_abs(ob)
ob.arg.is_abs = true;


% Cdiff1_power()
% for C.^2
function ob = Cdiff1_power(ob, p)
ob.arg.Cpower = ob.arg.Cpower * p;


% Cdiff1_ident_dup()
% y = I * x
%
function y = Cdiff1_ident_dup(arg, x)
y = x;


% Cdiff1_coef_powabs()
function coef = Cdiff1_coef_powabs(coef, arg)
coef = coef .^ arg.Cpower;
if arg.is_abs
	coef = abs(coef);
end


% Cdiff1_i1()
% indices for 1st-order differences
function ii = Cdiff1_i1(offset, nn)
ii = (1 + max(offset,0)) : (nn + min(offset,0));


% Cdiff1_i2()
% indices for 2nd-order differences
function ii = Cdiff1_i2(offset, nn)
offset = abs(offset);
ii = (1 + offset) : (nn - offset);


% Cdiff1_mex_forw()
% y = C * x
% both [(N) *L] (handles multi)
%
function y = Cdiff1_mex_forw(arg, x)

if arg.is_abs
	diff_str = sprintf('diff%d,forwA', arg.order);
else
	diff_str = sprintf('diff%d,forw%d', arg.order, arg.Cpower);
end

if isreal(x)
	y = penalty_mex(diff_str, single(x), arg.offset, arg.isize32);
else
	yr = penalty_mex(diff_str, single(real(x)), arg.offset, arg.isize32);
	yi = penalty_mex(diff_str, single(imag(x)), arg.offset, arg.isize32);
	y = complex(yr, yi);
end
y = reshapee(y, arg.isize, []); % [(N) 1 *L] to [(N) *L] due to mex


% Cdiff1_mex_forw_Fatrix()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_mex_forw_Fatrix(arg, x)

[x ei] = embed_in(x, true([arg.isize 1]), arg.nn); % [(N) *L]
y = Cdiff1_mex_forw(arg, x);
y = ei.shape(y); % [*N (L)] or [(N) (L)]


% Cdiff1_mex_back()
% x = C' * y
% both [(N) *L] (handles multi)
%
function x = Cdiff1_mex_back(arg, y)

if arg.is_abs
	diff_str = sprintf('diff%d,backA', arg.order);
else
	diff_str = sprintf('diff%d,back%d', arg.order, arg.Cpower);
end

y = reshapee(y, arg.isize, 1, []); % [(N) *L] to [(N) 1 *L] for mex
if isreal(y)
	x = penalty_mex(diff_str, single(y), arg.offset, arg.isize32);
else
	xr = penalty_mex(diff_str, single(real(y)), arg.offset, arg.isize32);
	xi = penalty_mex(diff_str, single(imag(y)), arg.offset, arg.isize32);
	x = complex(xr, xi);
end


% Cdiff1_mex_back_Fatrix()
% x = C' * y
% both [*N (L)] or [(N) (L)]
%
function x = Cdiff1_mex_back_Fatrix(arg, y)

[y eo] = embed_out(y, arg.isize); % [(N) *L]
x = Cdiff1_mex_back(arg, y);
x = eo.shape(x, true([arg.isize 1]), arg.nn); % [*N (L)] or [(N) (L)]


%% circshift

% Cdiff1_cs_forw()
% y = C * x
% both [(N) *L] (does_many)
%
function y = Cdiff1_cs_forw(arg, x)

displace = arg.displace;

switch arg.order
case 1
	coef = Cdiff1_coef_powabs(-1, arg);
	y = x + coef * circshift(x, displace);

case 2
	coef = Cdiff1_coef_powabs([2 -1], arg);
	y = coef(1) * x ...
		+ coef(2) * circshift(x, displace) ...
		+ coef(2) * circshift(x, -displace);

otherwise
	fail 'bug'
end


% Cdiff1_cs_forw_Fatrix()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_cs_forw_Fatrix(arg, x)

[x ei] = embed_in(x, true([arg.isize 1]), arg.nn); % [(N) *L]
y = Cdiff1_cs_forw(arg, x);
y = ei.shape(y); % [*N (L)] or [(N) (L)]


% Cdiff1_cs_back()
% x = C' * y
% both [(N) *L]
%
function x = Cdiff1_cs_back(arg, y)
arg.displace = -arg.displace; % trick
x = Cdiff1_cs_forw(arg, y);


% Cdiff1_cs_back_Fatrix()
% x = C' * y
% both [*N (L)] or [(N) (L)]
%
function x = Cdiff1_cs_back_Fatrix(arg, y)
arg.displace = -arg.displace; % trick
x = Cdiff1_cs_forw_Fatrix(arg, y);


%% diff


% Cdiff1_diff_forw()
% y = C * x
% both [(N) *L] (does_many)
%
function y = Cdiff1_diff_forw(arg, x)

if arg.order ~= 1 || arg.Cpower ~= 1 || arg.is_abs ~= 0
	fail 'only order=1 power=1 is_abs=0 done for "diff"'
end
y = diff(x, 1, arg.diff_dim);
% non-periodic end conditions:
siz = size(x);
siz(arg.diff_dim) = 1;
y = cat(arg.diff_dim, zeros(siz, class(y)), y);


% Cdiff1_diff_back()
% x = C' * y
% both [(N) *L]
%
function x = Cdiff1_diff_back(arg, y)
if arg.diff_dim ~= 1
	fail 'only 1st dim done'
end
switch numel(arg.isize)
case 1
	y([1 end+1],:) = 0;
case 2
	y([1 end+1],:,:) = 0;
case 3
	y([1 end+1],:,:,:) = 0;
otherwise
	fail 'only up to 3d done'
end

x = -diff(y, 1, arg.diff_dim);


%% for1


% Cdiff1_for1_loop1()
% a small function for matlab to compile "just in time"
% both [(N)]
function out = Cdiff1_for1_loop1(in, i1, i2, off, coef)
out = zeros(size(in), class(in));
switch class(in)
case 'single'
	coef = single(coef);
case 'double'
	coef = double(coef);
otherwise
	fail 'not done'
end
for ii = i1:i2
	out(ii) = in(ii) + coef * in(ii-off);
end


% Cdiff1_for1_loop2()
% a small function for matlab to compile "just in time"
% both [(N)]
function out = Cdiff1_for1_loop2(in, i1, i2, off, coef)
out = zeros(size(in), class(in));
switch class(in)
case 'single'
	coef = single(coef);
case 'double'
	coef = double(coef);
otherwise
	fail 'not done'
end
for ii = i1:i2
	out(ii) = coef(1) * in(ii) + coef(2) * (in(ii-off) + in(ii+off));
end


% Cdiff1_for1_forw()
% y = C * x
% both [(N)]
%
function y = Cdiff1_for1_forw(arg, x)

off = arg.offset;
switch arg.order
case 1
	coef = Cdiff1_coef_powabs(-1, arg);

	if 1 % faster!
		y = Cdiff1_for1_loop1(x, ...
			1+max(off,0), arg.nn+min(off,0), off, coef);
	else % slower!
		y = zeros(size(x), class(x));
		for ii = Cdiff1_i1(off, arg.nn)
			y(ii) = x(ii) + coef * x(ii-off);
		end
	end

case 2
	coef = Cdiff1_coef_powabs([2 -1], arg);
	if 1 % faster!
		y = Cdiff1_for1_loop2(x, ...
			1+abs(off), arg.nn-abs(off), off, coef);
	else % slower!
		y = zeros(size(x), class(x));
		for ii = Cdiff1_i2(off, arg.nn)
			y(ii) = coef(1) * x(ii) + coef(2) * (x(ii-off) + x(ii+off));
		end
	end

otherwise
	fail 'bug'
end


% Cdiff1_for1_forw_Fatrix()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_for1_forw_Fatrix(arg, x)

flag_array = 0;
if size(x,1) ~= arg.nn
	diml = size(x); diml = diml((1+length(arg.isize)):end);
	x = reshapee(x, arg.nn, []); % [*N *L]
	flag_array = 1;
end

y = zeros(size(x), class(x));
for ll = 1:ncol(x)
	y(:,ll) = Cdiff1_for1_forw(arg, x(:,ll));
end

if flag_array
	y = reshapee(y, arg.isize, diml); % [(N) (L)]
end


% Cdiff1_ind_forw()
% y = C * x
% both [(N)]
%
function y = Cdiff1_ind_forw(arg, x)

off = arg.offset;
y = zeros(size(x), class(x));

switch arg.order
case 1
	coef = Cdiff1_coef_powabs(-1, arg);
	ii = Cdiff1_i1(arg.offset, arg.nn);

	if 0 % timings: indexing slows it down a lot!
		cpu etic
		y = x + coef * x;
		cpu etoc 'with no index'

		cpu etic
		y = x(1:end-1) + x(2:end);
		cpu etoc 'with wired index'

		cpu etic
		y = x(ii) + x(ii - off);
		cpu etoc 'with index and no coef'

		cpu etic
		y(ii) = x(ii) + coef * x(ii-off);
		cpu etoc 'with index and coef'
	end

	y(ii) = x(ii) + coef * x(ii-off);

case 2
	coef = Cdiff1_coef_powabs([2 -1], arg);
	ii = Cdiff1_i2(arg.offset, arg.nn);
	y(ii) = coef(1) * x(ii) + coef(2) * (x(ii-off) + x(ii+off));

otherwise
	fail bug
end


% Cdiff1_ind_forw_Fatrix()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_ind_forw_Fatrix(arg, x)

flag_array = 0;
if size(x,1) ~= arg.nn
	diml = size(x); diml = diml((1+length(arg.isize)):end);
	x = reshapee(x, arg.nn, []); % [*N *L]
	flag_array = 1;
end

y = zeros(size(x), class(x));
for ll=1:ncol(x)
	y(:,ll) = Cdiff1_ind_forw(arg, x(:,ll));
end

if flag_array
	y = reshapee(y, arg.isize, diml); % [(N) (L)]
end


% Cdiff1_ind_back()
% x = C' * y
% both [(N)]
%
function x = Cdiff1_ind_back(arg, y)

off = arg.offset;
x = zeros(size(y), class(y));

switch arg.order
case 1
	coef = Cdiff1_coef_powabs(-1, arg);
	ii = Cdiff1_i1(arg.offset, arg.nn);
%	y(ii) = coef(1) * x(ii) + coef(2) * x(ii - arg.offset);
	x(ii) = y(ii);
	x(ii-off) = x(ii-off) + coef * y(ii);

case 2
	coef = Cdiff1_coef_powabs([2 -1], arg);
	ii = Cdiff1_i2(arg.offset, arg.nn);
%	y(ii) = coef(1) * x(ii) + coef(2) * (x(ii-off) + x(ii+off));
	x(ii) = coef(1) * y(ii);
	x(ii-off) = x(ii-off) + coef(2) * y(ii);
	x(ii+off) = x(ii+off) + coef(2) * y(ii);
otherwise
	fail bug
end


% Cdiff1_ind_back_Fatrix()
% x = C' * y
% both [*N (L)] or [(N) (L)]
%
function x = Cdiff1_ind_back_Fatrix(arg, y)

flag_array = 0;
if size(y,1) ~= arg.nn
	diml = size(y); diml = diml((1+length(arg.isize)):end);
	y = reshapee(y, arg.nn, []); % [*N *L]
	flag_array = 1;
end

x = zeros(size(y), class(y));
for ll=1:ncol(y)
	x(:,ll) = Cdiff1_ind_back(arg, y(:,ll));
end

if flag_array
	x = reshapee(x, arg.isize, diml); % [(N) (L)]
end


% Cdiff1_coef_setup()
% precompute filter coefficients for 'circshfit' 'convn' 'imfilter' versions
%
function arg = Cdiff1_coef_setup(arg)

dd = arg.displace;
cdim = 1 + 2 * abs(dd); % trick: odd always to help adjoint and edges
coef = zeros([cdim 1]);

switch arg.order
case 0
	coef = 1;

case 1 % [0 1 -1]
	mid = 1 + abs(dd);
	tmp = num2cell(mid);
	coef(tmp{:}) = 1;
	tmp = num2cell(mid + dd);
	coef(tmp{:}) = -1;

case 2 % [-1 2 -1]
	mid = 1 + abs(dd);
	tmp = num2cell(mid);
	coef(tmp{:}) = 2;
	tmp = num2cell(mid + dd);
	coef(tmp{:}) = -1;
	tmp = num2cell(mid - dd);
	coef(tmp{:}) = -1;

otherwise
	fail('order')
end

arg.coef = coef;


% Cdiff1_filt_setup()
% precompute filter coefficients for 'convn' and 'imfilter' versions
%
function arg = Cdiff1_filt_setup(arg)

switch arg.type_diff
case 'convn'
	arg.filt_forw = @(x, c) convn(x, c, 'same');
%		padn(convn(x, c, 'valid'), arg.size);
case 'imfilter'
	arg.filt_forw = @(x, c) imfilter(x, c, 'circular', 'conv', 'same');
otherwise
	fail 'bug'
end


% Cdiff1_filt_forw_Fatrix()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_filt_forw_Fatrix(arg, x)

nd = length(arg.isize);

[x ei] = embed_in(x, true([arg.isize 1]), arg.nn); % [(N) *L]

coef = Cdiff1_coef_powabs(arg.coef, arg);

LL = size(x, numel(arg.isize)+1);
if LL == 1
	y = arg.filt_forw(x, coef);
else
	y = zeros(arg.nn, LL, class(x));
	for ll=1:LL
		tmp = stackpick(x,ll);
		tmp = arg.filt_forw(tmp, coef);
		y(:,ll) = tmp(:);
	end
	y = reshape(y, [arg.isize LL]);
end

y = ei.shape(y); % [*N (L)] or [(N) (L)]


% Cdiff1_filt_back_Fatrix()
% x = C' * y
% both [*N (L)] or [(N) (L)]
%
function x = Cdiff1_filt_back_Fatrix(arg, y)
arg.coef = flipdims(arg.coef); % trick
x = Cdiff1_filt_forw_Fatrix(arg, y);


% Cdiff1_sp_setup()
% setup for sparse version
%
function C = Cdiff1_sp_setup(arg)

if arg.order == 1
	off = arg.offset;
	ii = Cdiff1_i1(arg.offset, arg.nn);
	C	= sparse(ii, ii, 1, arg.nn, arg.nn) ...
		- sparse(ii, ii-arg.offset, 1, arg.nn, arg.nn);

else % order = 2
	ii = Cdiff1_i2(arg.offset, arg.nn);
	C	= sparse(ii, ii, 2, arg.nn, arg.nn) ...
		- sparse(ii, ii-arg.offset, 1, arg.nn, arg.nn) ...
		- sparse(ii, ii+arg.offset, 1, arg.nn, arg.nn);
end


% Cdiff1_sp_forw()
% y = C * x
% both [*N (L)] or [(N) (L)]
%
function y = Cdiff1_sp_forw(arg, x)

C = arg.C;
if arg.is_abs
	C = abs(C);
end
if arg.Cpower ~= 1
	C = C .^ arg.Cpower;
end

y = C * x;


% Cdiff1_sp_back()
% x = C' * y
% both [*N (L)] or [(N) (L)]
%
function x = Cdiff1_sp_back(arg, y)

C = arg.C;
if arg.is_abs
	C = abs(C);
end
if arg.Cpower ~= 1
	C = C .^ arg.Cpower;
end

x = C' * y;
