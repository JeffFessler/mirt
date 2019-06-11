 function ob = Gtomo2_wtmex(arg1, varargin)
%function ob = Gtomo2_wtmex(sg, ig, options)
%function ob = Gtomo2_wtmex(file, options)
%function ob = Gtomo2_wtmex(arg_pairs, options)
%|
%| Generate a 2D tomographic system model based on aspire's wtfmex() routine.
%|
%| in
%|	sg	strum		sino_geom()
%|	ig	strum		image_geom()
%|
%| options
%|	'grouped'		row or col (default: row)
%|	'nthread'		pthreads for multiple processors (default: 1)
%|	'mask'	[nx ny]		logical support mask, has precedence
%|	'pairs' {}		cell array of name/value pairs for aspire_pair()
%|	'class'			'Fatrix' (default) or 'fatrix2'
%|
%| out
%|	ob	[nb*na np]	Fatrix system "matrix", np = sum(ig.mask(:))
%|
%| Copyright 2006-3-3, Jeff Fessler, University of Michigan

if nargin == 1 && streq(arg1, 'test'), Gtomo2_wtmex_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% options
arg.grouped = 'row';
arg.nthread = 1;
arg.mask = [];
arg.chat = 0;
arg.class = 'Fatrix'; % todo: keep until block version fixed!
%arg.class = 'fatrix2';
arg.pairs = {};


if ischar(arg1) % given wtf or argument pairs

	arg = vararg_pair(arg, varargin);

	if exist(arg1, 'file') % file.wtf
		if ~isempty(arg.mask), error 'no mask option given .wtf', end
		arg.file = arg1;
		[arg.buff nx ny nb na] = ...
			wtfmex('asp:read', arg.file, int32(arg.chat));

	else % arg_pair
		arg.arg = arg1;
		if isempty(arg.mask)
			mask_arg = {};
		else
			if ~islogical(arg.mask), error 'need logical mask', end
			mask_arg = uint8(arg.mask);
		end

		arg.buff = wtfmex('asp:gensys', arg.arg', arg.grouped, mask_arg{:});

		nx = str2num(arg_get(arg.arg, 'nx'));
		ny = str2num(arg_get(arg.arg, 'ny'));
		nb = str2num(arg_get(arg.arg, 'nb'));
		na = str2num(arg_get(arg.arg, 'na'));
	end

	mask = wtfmex('asp:mask', arg.buff) > 0;
	arg.ig = image_geom('nx', nx, 'ny', ny, 'dx', -1, 'mask', mask);
	arg.sg = sino_geom('par', 'nb', nb, 'na', na, ...
		'strip_width', 'dr'); % unused but prevents warning

else % given sino_geom() and image_geom()

	arg.sg = arg1;
	arg.ig = varargin{1};
	varargin = {varargin{2:end}};
	arg = vararg_pair(arg, varargin);
	if ~isempty(arg.mask)
		arg.ig.mask = arg.mask;
	end

	arg.aspire_arg = aspire_pair(arg.sg, arg.ig, 'support', 'array', ...
		arg.pairs{:}); % trick: pairs can override default support

	arg.buff = wtfmex('asp:gensys', arg.aspire_arg', arg.grouped, uint8(arg.ig.mask));

end

arg.power = 1;
arg.np = sum(arg.ig.mask(:));

arg.nthread = int32(arg.nthread);

% wtfmex() method(s)
arg.stayman2_factors = @(G, wi) wtfmex('asp:stayman2', G.arg.buff, single(wi));
arg.nuyts2_factors = @(G, wi) wtfmex('asp:nuyts2', G.arg.buff, single(wi));

arg.odim = [arg.sg.nb arg.sg.na];

switch arg.class

case 'fatrix2'
	forw = @(arg, x) ...
		wtfmex('asp:forw', arg.buff, arg.nthread, single(x), int32(arg.chat));
	back = @(arg, y) Gtomo2_wtmex_back_fatrix2(arg, y);

	forw_block = @(arg, x, iblock, nblock) ...
		Gtomo2_wtmex_block_forw(arg, [], x, iblock, nblock);

	back_block = @(arg, y, iblock, nblock) ...
		Gtomo2_wtmex_block_back(arg, [], y, iblock, nblock);

	ob = fatrix2('arg', arg, ...
		'idim', arg.ig.dim, 'odim', arg.odim, ...
		'imask', arg.ig.mask, ...
		'abs', @(ob) ob, ... % nonnegative
		'forw', forw, 'back', back, ...
		'forw_block', forw_block, 'back_block', back_block, ...
		'free', @Gtomo2_wtmex_free, ...
		'power', @Gtomo2_wtmex_power);

case 'Fatrix'
	arg.nd = arg.sg.nb * arg.sg.na;
	dim = [arg.nd arg.np]; % trick: make it masked by default!
	ob = Fatrix(dim, arg, 'caller', 'Gtomo2_wtmex', ...
		'abs', @(ob) ob, ... % nonnegative
		'forw', @Gtomo2_wtmex_forw, 'back', @Gtomo2_wtmex_back, ...
		'free', @Gtomo2_wtmex_free, ...
		'mtimes_block', @Gtomo2_wtmex_mtimes_block, ...
		'power', @Gtomo2_wtmex_power);

otherwise
	fail('unknown class "%s"', arg.class)
end

if arg.chat
	wtfmex('asp:print', arg.buff);
end


% Gtomo2_wtmex_back_fatrix2(): x = G' * y
function x = Gtomo2_wtmex_back_fatrix2(arg, y)
if arg.power == 2
	x = wtfmex('asp:back2', arg.buff, arg.nthread, single(y), int32(arg.chat));
elseif arg.power == 1
	x = wtfmex('asp:back', arg.buff, arg.nthread, single(y), int32(arg.chat));
else
	fail('power=%d not done', arg.power)
end


% Gtomo2_wtmex_forw(): y = G * x
function y = Gtomo2_wtmex_forw(arg, x)
y = Gtomo2_wtmex_mtimes_block(arg, 0, x, 1, 1);


% Gtomo2_wtmex_back(): x = G' * y
function x = Gtomo2_wtmex_back(arg, y)
x = Gtomo2_wtmex_mtimes_block(arg, 1, y, 1, 1);


% Gtomo2_wtmex_mtimes_block()
function y = Gtomo2_wtmex_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	fun = @Gtomo2_wtmex_block_back;
else
	fun = @Gtomo2_wtmex_block_forw;
end
y = embed_mult(fun, arg, is_transpose, x, istart, nblock, ...
	arg.ig.mask, arg.np, [arg.sg.nb arg.sg.na], 1);


% Gtomo2_wtmex_block_forw()
function y = Gtomo2_wtmex_block_forw(arg, dummy, x, istart, nblock)

if arg.power ~= 1, fail('power=%d not done', arg.power), end
% todo: could square the values in arg.buff !

if nblock == 1
	y = wtfmex('asp:forw', arg.buff, arg.nthread, single(x), int32(arg.chat));
else
	y = wtfmex('asp:proj,block', arg.buff, arg.nthread, single(x), ...
			int32(istart-1), int32(nblock), int32(arg.chat));

	% todo: add warning!
	% fix: todo extract the relevant columns - should do in wtfmex?
	ia = istart:nblock:arg.sg.na;
	y = y(:,ia,:);
end


% Gtomo2_wtmex_block_back()
function x = Gtomo2_wtmex_block_back(arg, dummy, y, istart, nblock)

if nblock == 1
	if arg.power == 2
		x = wtfmex('asp:back2', arg.buff, arg.nthread, single(y), int32(arg.chat));
	elseif arg.power == 1
		x = wtfmex('asp:back', arg.buff, arg.nthread, single(y), int32(arg.chat));
	else
		fail('power=%d not done', arg.power)
	end
else
	if arg.power ~= 1, fail('power=%d not done', arg.power), end
	% the following block might be needed for fatrix2, but NOT for Fatrix
	% because for Fatrix this is built into embed_mult.m
	if 0 && nblock ~= 1 % todo: build into wtfmex
		tmp = zeros(arg.odim);
		ia = istart:nblock:arg.sg.na;
		tmp(:,ia,:) = y;
		y = tmp;
	end
	x = wtfmex('asp:back,block', arg.buff, arg.nthread, single(y), ...
			int32(istart-1), int32(nblock), int32(arg.chat));
end


% Gtomo2_wtmex_free()
function Gtomo2_wtmex_free(arg)
printm 'freeing wtfmex static memory'
%% wtfmex('free');


% Gtomo2_wtmex_power()
function ob = Gtomo2_wtmex_power(ob, sup)
ob.arg.power = ob.arg.power * sup;
