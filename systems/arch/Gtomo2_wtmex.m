 function ob = Gtomo2_wtmex(arg1, varargin)
%function ob = Gtomo2_wtmex(sg, ig, options)
%function ob = Gtomo2_wtmex(file, options)
%function ob = Gtomo2_wtmex(arg_pairs, options)
%
% Generate a tomographic system model based on aspire's wtfmex() routine.
%
% in
%	sg	strum		sino_geom()
%	ig	strum		image_geom()
%
% options
%	'grouped'		row or col (default: row)
%	'nthread'		pthreads for multiple processors (default: 1)
%	'mask'	[nx,ny]		logical support mask, has precedence
%	'pairs' {}		cell array of name/value pairs for aspire_pair()
%
% out
%	ob	[nb*na,np]	Fatrix system "matrix", np = sum(ig.mask(:))
%
% Copyright 2006-3-3, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(arg1, 'test'), Gtomo2_wtmex_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% options
arg.grouped = 'row';
arg.nthread = 1;
arg.mask = [];
arg.chat = 0;
arg.pairs = {};

% given wtf or argument pairs
if ischar(arg1)

	arg = vararg_pair(arg, varargin);

	if exist(arg1, 'file') % file.wtf
		if ~isempty(arg.mask), error 'no mask option given .wtf', end
		arg.file = arg1;
		[nx ny nb na] = wtfmex('read', arg.file);

	else % arg_pair
		arg.arg = arg1;
		if isempty(arg.mask)
			mask_arg = {};
		else
			if ~islogical(arg.mask), error 'need logical mask', end
			mask_arg = uint8(arg.mask);
		end

		wtfmex('gensys', arg.arg', arg.grouped, mask_arg{:});

		nx = str2num(arg_get(arg.arg, 'nx'));
		ny = str2num(arg_get(arg.arg, 'ny'));
		nb = str2num(arg_get(arg.arg, 'nb'));
		na = str2num(arg_get(arg.arg, 'na'));
	end

	arg.ig = image_geom('nx', nx, 'ny', ny, 'dx', -1, ...
		'mask', wtfmex('mask') > 0);
	arg.sg = sino_geom('par', 'nb', nb, 'na', na);

% given sino_geom() and image_geom()
else

	arg.sg = arg1;
	arg.ig = varargin{1};
	varargin = {varargin{2:end}};
	arg = vararg_pair(arg, varargin);
	if ~isempty(arg.mask)
		arg.ig.mask = arg.mask;
	end

	arg.aspire_arg = aspire_pair(arg.sg, arg.ig, 'support', 'array', ...
		arg.pairs{:}); % trick: pairs can override default support

	wtfmex('gensys', arg.aspire_arg', arg.grouped, uint8(arg.ig.mask));

end

arg.power = 1;
arg.nd = arg.sg.nb * arg.sg.na;
arg.np = sum(arg.ig.mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtomo2_wtmex', ...
	'forw', @Gtomo2_wtmex_forw, 'back', @Gtomo2_wtmex_back, ...
	'free', @Gtomo2_wtmex_free, ...
	'mtimes_block', @Gtomo2_wtmex_mtimes_block);

if arg.chat
	wtfmex('print');
end


%
% Gtomo2_wtmex_forw(): y = G * x
%
function y = Gtomo2_wtmex_forw(arg, x)
y = Gtomo2_wtmex_mtimes_block(arg, 0, x, 1, 1);


%
% Gtomo2_wtmex_back(): x = G' * y
% full backprojection
%
function x = Gtomo2_wtmex_back(arg, y)
if arg.power == 2
	[flag_column y] = Gtomo2_wtmex_y_shape(arg, y, istart, nblock);
	x = wtfmex('chat', arg.chat, 'back2', single(y));
	if flag_column, x = arg.ig.maskit(x); end
else
	x = Gtomo2_wtmex_mtimes_block(arg, 1, y, 1, 1);
end


%
% Gtomo2_wtmex_mtimes_block()
%
function y = Gtomo2_wtmex_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	y = Gtomo2_wtmex_block_back(arg, x, istart, nblock);
else
	y = Gtomo2_wtmex_block_forw(arg, x, istart, nblock);
end


%
% Gtomo2_wtmex_block_forw()
%
function y = Gtomo2_wtmex_block_forw(arg, x, istart, nblock)

if arg.power ~= 1, error('power=%d not done', arg.power), end

% if needed, convert concise column to 3d array
flag_column = 0;
if size(x,1) == arg.np
	flag_column = 1;
	x = arg.ig.embed(x);
end

if nblock == 1
	y = wtfmex('chat', arg.chat, 'mult', single(x));
else
	y = wtfmex('chat', arg.chat, 'proj,block', single(x), ...
			int32(istart-1), int32(nblock));

	% fix: extract the relevant columns - should do in wtfmex?
	ia = istart:nblock:arg.sg.na;
	y = y(:,ia,:);
end

y = double6(y);

if flag_column % column in yields column out.
	y = reshape(y, arg.sg.nb*size(y,2), []);
end


%
% Gtomo2_wtmex_block_back()
%
function x = Gtomo2_wtmex_block_back(arg, y, istart, nblock)

if arg.power ~= 1, error('power=%d not done', arg.power), end

[flag_column y] = Gtomo2_wtmex_y_shape(arg, y, istart, nblock);
if nblock == 1
	x = wtfmex('chat', arg.chat, 'back', single(y));
else
	x = wtfmex('chat', arg.chat, 'back,block', single(y), ...
			int32(istart-1), int32(nblock));
end
x = double6(x);

if flag_column, x = arg.ig.maskit(x); end


%
% Gtomo2_wtmex_y_shape(arg, y, istart, nblock)
%
function [flag_column, y] = Gtomo2_wtmex_y_shape(arg, y, istart, nblock)

flag_column = 0;
nb = arg.sg.nb;
na = arg.sg.na;

if nblock == 1
	if size(y,1) == nb * na % [nb*na,?]
		flag_column = 1;
		y = arg.sg.shape(y); % [nb,na,?]
	end
else
	ia = istart:nblock:na;
	nv = length(ia);
	if size(y,1) == nb * nv % [nb*nv,?]
		flag_column = 1;
		tmp = zeros(nb, na, size(y,2));
		tmp(:,ia,:) = reshape(y, nb, nv, []);
		y = reshape(tmp, nb*na, []); % [nb*na,?]
	elseif size(y,1) == nb
		error todo
	else
		error bug
	end
end


%
% Gtomo2_wtmex_free()
%
function Gtomo2_wtmex_free(arg)
printm('freeing wtfmex static memory')
wtfmex('free');



%
% Gtomo2_wtmex_power()
%
function ob = Gtomo2_wtmex_power(ob, sup)
ob.arg.power = ob.arg.power * sup;
