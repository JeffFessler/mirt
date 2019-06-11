 function ob = Gtomo3(sys_type, mask, nx, ny, nz, varargin)
%function ob = Gtomo3(sys_type, mask, nx, ny, nz, options)
%|
%| Construct Gtomo3 object, for 3D tomographic image reconstruction,
%| especially for SPECT and perhaps PET.  (For CT use Gcone; see below.)
%| This object does 3d forward and backprojection using mex files.
%| See Gtomo3_test.m for example usage.
%|
%| Basically, you create an object calling:
%|		A = Gtomo3(sys_type, mask)
%| and then you can use it thereafter by typing commands like
%|		y = A * x;
%| which will auto-magically do the multiplication by calling a mex file.
%|
%| in
%|	sys_type	string	see 3D ASPIRE User's Guide
%|	mask	[nx ny nz] or string
%|		mask can be a filename, in which case the file is read.
%|		mask can be empty (for some system types, including 3s),
%|		in which case it is determined from backprojecting a "1" vector.
%|		Else, mask must be [nx ny nz] binary array of object support.
%|
%| options
%|	'nthread' 1, 2...	1, unless you have a multi-processor system
%|	'chat' 0|1		verbosity.  default: 1.
%|	'checkmask' 0|1		verify that the mask is correct.  default: 1.
%|	'view2d' 0|1		for 2z and 2dsc system types, use 2d views
%|				instead of sinograms, for use with OS.
%|				default: 0
%|	'permute213' 0|1	convert from 123 to 213 output order
%|				(e.g., to overcome vu vs st order of 3l@)
%|				default: 0
%|	'class' ''		default 'fatrix2'
%|	'scale' [1]		scale factor.  default '' (none)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|
%| For more help, type 'f3d_mex' and see Gtomo3_test.m
%|
%| For X-ray CT, using Gcone() is recommended instead of Gtomo3().
%| Nevertheless, the following alternate usage is permitted for CT here:
%| function ob = Gtomo3(cg, ig, options)
%|
%| Copyright 01-04-22, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sys_type, 'test'), Gtomo3_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% handle new usage (for CT)
% (cg, ig, [options])
if isa(sys_type, 'strum')
	cg = sys_type;
	ig = mask;
	options = {};
	if isvar('nx')
		options = {nx};
		if isvar('ny')
			options = {nx, ny};
			if isvar('nz')
				options = {nx, ny, nz};
			end
		end
	end
	options = {options{:}, varargin{:}};

	sys_type = aspire_pair(cg, ig, 'system', '3l');
	ob = Gtomo3(sys_type, ig.mask, ig.nx, ig.ny, ig.nz, options{:});
return
end

arg.sys_type = sys_type;
arg.nx = nx;
arg.ny = ny;
arg.nz = nz;
arg.nxyz = nx*ny*nz;
arg.mask = logical(mask);

% defaults
arg.nthread = 1;
arg.chat = true;
arg.checkmask = false;
arg.scale = [];
arg.view2d = false;
arg.permute213 = false;
arg.class = 'fatrix2';
% arg.class = 'fatrix2';
arg.nthread_max = 8; % > 8 threads may crash matlab in Gtomo3_test for '3s'
		% due to f3d_mex bugs?  increase this at your own risk!

if ischar(arg.mask)
	arg.mask = fld_read(arg.mask); % read mask file
end

% options
arg = vararg_pair(arg, varargin);

if arg.nthread > arg.nthread_max
	warn('%d threads downgraded to %d, proceed with caution!', ...
		arg.nthread, arg.nthread_max)
	arg.nthread = arg.nthread_max;
end

if arg.view2d && ~streq(arg.sys_type, '2', 1)
	fail 'view2d option only relevant to 2dsc and 2z system types'
end
if arg.permute213 && ~streq(arg.sys_type, '3l', 2)
	warn 'permute213 designed only for "3l" system model!?'
end

% if no mask, then backproject to find it!
if isempty(arg.mask)
	arg.optr = f3d_mex('init', arg.sys_type, uint8(ones(nx,ny,nz)), ...
		int32(arg.nx), int32(arg.ny), int32(arg.nz), ...
		int32(arg.nthread), int32(arg.chat))
	warn('creating mask by summing will be a bit slow')
	t = double(f3d_mex('dims', arg.optr)); % data dimensions, matlab wants doubles
	arg.mask = f3d_mex('back', ones(t, 'single'), ...
			int32(0), int32(1), arg.optr, int32(0)) > 0;
	f3d_mex('free:1', arg.optr)
end

% initialize mex file
arg.optr = f3d_mex('init', arg.sys_type, uint8(arg.mask), ...
	int32(arg.nx), int32(arg.ny), int32(arg.nz), ...
	int32(arg.nthread), int32(arg.chat));
if arg.chat
%	printm('arg.optr = %lX', arg.optr)
end

arg.odim = Gtomo3_nn(arg);
arg.nd = prod(arg.odim);
arg.np = sum(arg.mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

% trick: handle 3l@ uv vs st mismatch (if requested) by permuting.
if arg.permute213
%	if ~isempty(after), fail('permute213 and scale not done'), end
%	after = @(x, is_transpose, istart, nblock) permute(x, [2 1 3]);
%	after = {'cascade_after', after};
	warn 'still testing permute213'
end

abs_arg = {'abs', @(ob) ob}; % aij nonnegative, except perhaps SPECT rotate
if streq(arg.sys_type, '3s@', 3)
	abs_arg = {}; % todo: check for spline_filter 3, else ok
end

switch arg.class
case 'fatrix2'
	forw = @(arg, x) Gtomo3_block_forw_array(arg, x, 1, 1);
	back = @(arg, y) Gtomo3_block_back_array(arg, y, 1, 1);
	forw_block = @(arg, x, iblock, nblock) ...
		Gtomo3_block_forw_array(arg, x, iblock, nblock);
	back_block = @(arg, y, iblock, nblock) ...
		Gtomo3_block_back_array(arg, y, iblock, nblock);

	ob = fatrix2('mask', arg.mask, 'arg', arg, abs_arg{:}, ...
		'does_many', 1, ...
		'forw', forw, 'back', back, ...
		'forw_block', forw_block, 'back_block', back_block);
	if ~isempty(arg.scale)
		ob = arg.scale * ob;
		arg.scale = [];
	end

case 'Fatrix' % build Fatrix object
	if ~isempty(arg.scale)
		after = {'cascade_after', arg.scale};
	else
		after = {};
	end

	ob = Fatrix(dim, arg, 'caller', 'Gtomo3', abs_arg{:}, ...
		'forw', @Gtomo3_forw_Fatrix, ...
		'back', @Gtomo3_back_Fatrix, ...
		'free', @Gtomo3_free, ...
		'mtimes_block', @Gtomo3_mtimes_block_Fatrix, ...
		after{:});

otherwise
	fail('class %s', arg.class)
end

if arg.chat
	f3d_mex('show:1', arg.optr, int32(arg.chat)); % display header info
end

% check mask
if arg.checkmask
	printm 'Warn checking mask, which may be unnecessary'
	tmp = sum(ob).' > 0;
	tmp = embed(tmp, arg.mask);
	if ~isequal(tmp, arg.mask)
		if exist('imfill') == 2
			tmpfill = imfill(tmp);
			if ~isequal(tmpfill, arg.mask)
				if im
					im plc 1 2
					im(1, arg.mask, 'input mask')
					im(2, tmp, 'computed mask')
				end
				warn('mask-sum mismatch not corrected by fill')
			else
				warn('mask-sum mismatch; corrected by filling')
			end
			printm('Warn: this geometry may have poor sampling')
		else
			fail 'mask inconsistent with sum'
		end
	end
	printm('mask checking completed')
end


% Gtomo3_nn(): data dimensions
function nn = Gtomo3_nn(arg)
nn = double(f3d_mex('dims', arg.optr)); % silly matlab wants size to be double
if arg.view2d
	nn = nn([1 3 2]); % swap projection angle and slice index
end
if arg.permute213
	nn = nn([2 1 3]); % trick: swap vu to st for 3l@
end


% Gtomo3_forw_Fatrix(): y = G * x
% full projection
function y = Gtomo3_forw_Fatrix(arg, x)
y = Gtomo3_mtimes_block_Fatrix(arg, 0, x, 1, 1);


% Gtomo3_back_Fatrix(): x = G' * y
% full backprojection
function x = Gtomo3_back_Fatrix(arg, y)
x = Gtomo3_mtimes_block_Fatrix(arg, 1, y, 1, 1);


% Gtomo3_mtimes_block_Fatrix()
function y = Gtomo3_mtimes_block_Fatrix(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = Gtomo3_block_back_Fatrix(arg, x, istart, nblock);
else
	y = Gtomo3_block_forw_Fatrix(arg, x, istart, nblock);
end


% Gtomo3_block_forw_Fatrix()
function y = Gtomo3_block_forw_Fatrix(arg, x, istart, nblock)

[x ei] = embed_in(x, arg.mask, arg.np);
y = Gtomo3_block_forw_array(arg, x, istart, nblock);
y = ei.shape(y);


% Gtomo3_block_back_Fatrix()
function x = Gtomo3_block_back_Fatrix(arg, y, istart, nblock)

na = arg.odim(end); % last index must be the one over which we subsetize
ia = istart:nblock:na;
[y eo] = embed_out(y, [arg.odim(1:end-1) length(ia)]); % [(M) *L]
x = Gtomo3_block_back_array(arg, y, istart, nblock);
x = eo.shape(x, arg.mask, arg.np);


% Gtomo3_free()
function Gtomo3_free(arg)
if arg.chat
	printm('freeing Gtomo3 object static memory')
end
f3d_mex('free:1', arg.optr);


%
% fatrix2 versions "_array"
%


% Gtomo3_block_forw_array()
% x: [nx ny nz *L]
% y: [odim(1:end-1) "na" *L] where "na" is # of angles in this block
function y = Gtomo3_block_forw_array(arg, x, istart, nblock)

for i4=1:size(x,4)
	y(:,:,:,i4) = f3d_mex('proj', single(x(:,:,:,i4)), ...
		int32(istart-1), int32(nblock), arg.optr, int32(arg.chat));
end
if arg.view2d % trick: for 2z and 2dsc
	y = permute(y, [1 3 2 4]);
end
if arg.permute213 % trick: for 3l@
	y = permute(y, [2 1 3 4]);
end
na = arg.odim(end); % last index must be the one over which we subsetize
ia = istart:nblock:na;
y = reshapee(y, prod(arg.odim(1:end-1)), na, []);
LL = size(y,3);
y = y(:,ia,:);
y = reshape(y, [arg.odim(1:end-1) length(ia) LL]); % [(M) *L]


% Gtomo3_block_back_array()
function x = Gtomo3_block_back_array(arg, y, istart, nblock)

na = arg.odim(end); % last index must be the one over which we subsetize
ia = istart:nblock:na;

n1 = prod(arg.odim(1:end-1));
y = reshapee(y, n1, length(ia), []);
yy = size(y); yy(2) = na; yy = zeros(yy);
yy(:,ia,:) = y;
dim = num2cell(arg.odim);
y = reshapee(yy, dim{:}, []);

if arg.permute213 % trick: for 3l@
	y = permute(y, [2 1 3 4]);
end
if arg.view2d % trick for 2d
	y = permute(y, [1 3 2 4]);
end
for i4 = 1:size(y,4)
	x(:,:,:,i4) = f3d_mex('back', single(y(:,:,:,i4)), ...
		int32(istart-1), int32(nblock), arg.optr, int32(arg.chat));
end
