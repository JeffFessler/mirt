 function ob = fatrix2_kroni(Mkron, ob)
%function ob = fatrix2_kroni(Mkron, ob)
%|
%| Construct fatrix2_block object of form kron(eye(Mkron), ob)
%| i.e., a kronecker product of identity matrix with object.
%|
%| in
%|	Mkron		natural number
%|	ob		fatrix2
%|
%| out
%|	ob		fatrix2 equivalent to kron(eye(Mkron), ob)
%|
%| Copyright 2011-09-29, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

arg.Mkron = Mkron;
arg.ob = ob;

arg.idim = [ob.idim arg.Mkron];
arg.odim = [ob.odim arg.Mkron];

if isempty(ob.imask)
	imask = [];
else
	imask = repmat(ob.imask, [ones(1,numel(ob.idim)) arg.Mkron]);
end

if isempty(ob.omask)
	omask = [];
else
	omask = repmat(ob.omask, [ones(1,numel(ob.odim)) arg.Mkron]);
end

if ~isempty(ob.handle_abs)
	abs_arg = {'abs', @fatrix2_kroni_abs};
else
	abs_arg = {};
end

% build object
tmp = sprintf('fatrix2_kronI, %s)', ob.caller);
ob = fatrix2('arg', arg, 'caller', tmp, ...
	'idim', arg.idim, 'odim', arg.odim, ...
	'imask', imask, 'omask', omask, ...
	abs_arg{:}, ...
	'forw', @fatrix2_kroni_forw, ...
	'back', @fatrix2_kroni_back);
%	'fatrix2_block', @fatrix2_kroni_block_setup, ...
%	'mtimes_block', @fatrix2_kroni_mtimes_block, ...


% fatrix2_kroni_forw(): y = A * x
function y = fatrix2_kroni_forw(arg, x)

y = cell(arg.Mkron, 1);
for ii=1:arg.Mkron
	tmp = stackpick(x,ii);
	y{ii} = arg.ob * tmp;
end
y = cat(numel(arg.odim), y{:});


% fatrix2_kroni_back(): x = A' * y
function x = fatrix2_kroni_back(arg, y)

o1 = arg.ob;
x = cell(arg.Mkron, 1);
for ii=1:arg.Mkron
	tmp = stackpick(y,ii);
	x{ii} = o1' * tmp;
	if numel(o1.odim) == 1 % trick for 1D odim cases
		x{ii} = fatrix2_embed(o1.imask, o1.idim, x{ii});
	end
end
x = cat(numel(arg.idim), x{:});


% fatrix2_kroni_abs()
function ob = fatrix2_kroni_abs(ob)
o1 = ob.arg.ob;
ob.arg.ob = abs(o1);


% fatrix2_kroni_block_setup()
% apply Gblock to the base block, to prepare it for later
function ob = fatrix2_kroni_block_setup(ob, varargin)
ob.arg.blocks = {Gblock(ob.arg.blocks{1}, ob.nblock)};
% todo: probably need to set up some other internal variables here
% for use inside fatrix2_kroni_mtimes_block()


% fatrix2_kroni_mtimes_block(): y = A{ib} * x
function y = fatrix2_kroni_mtimes_block(arg, is_transpose, x, iblock, nblock)

fail 'not done'

bl1 = arg.blocks{1}; % base block, already put through Gblock
warn 'todo: size check not done'

if ~is_transpose % forw
%	if nrow(x) ~= arg.dim(2)
%		fail('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
%	end
	y = [];
	for ii=1:arg.Mkron
		t = bl1{iblock} * x([arg.jstart(ii):arg.jend(ii)], :);
		y = [y; t];
	end

else % back
	y = x;
%	if nrow(y) ~= arg.dim(1), error 'bad y size', end
	x = [];
	for ii=1:arg.Mkron
		tmp = [arg.istart(ii):arg.iend(ii)]; % i list (if all data)
		% todo: we need a certain subset of that list
		fail('todo: kron subset backprojector not done')
		t = bl1{iblock}' * y(tmp, :);
		x = [x; t];
	end
	y = x;
end
