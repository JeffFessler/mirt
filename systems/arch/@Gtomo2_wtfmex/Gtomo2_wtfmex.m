 function ob = Gtomo2_wtfmex(file, chat, sp_group_type)
%function ob = Gtomo2_wtfmex(file, chat, sp_group_type)
%
% THIS IS OBSOLETE; USE Gtomo2_wtmex INSTEAD!
%
% Construct Gtomo2_wtfmex object, which does 2d forward and backprojection
% using the wtfmex cabilities for any system geometry.
% See Gtomo2_wtfmex_test.m for example usage.
% This object avoids the memory pigging property of matlab's sparse matrices
% using wtfmex's static system matrix capability.
% Basically, you create a wtfmex system matrix by calling:
%	G = Gtomo2_wtfmex(file)
% and then you can use it thereafter by typing commands like
%	y = G * x;
% which will auto-magically call wtfmex to do the multiplication
% in the compiled C mex program.
%
% Alternative usage: let "file" be geometry description from arg_pair().
%
% Copyright 2002-2-20, Jeff Fessler, The University of Michigan

if ~isvar('chat') | isempty(chat)
	chat = 0;
end

if ~isvar('sp_group_type') | isempty(sp_group_type)
	sp_group_type = 'row';
end

%
% create default object, as required by Mathworks
%
ob = Gtomo2(chat);
ob.file = '';
ob.arg = '';

if ~nargin % required by Mathworks
	ob = class(ob, 'Gtomo2_wtfmex');
return
end

warning 'Gtomo2_wtfmex is obsolete.  use Gtomo2_wtmex instead'

if exist(file, 'file')
	ob.file = file;

	[ob.nx, ob.ny, ob.nb, ob.na] = wtfmex('read', ob.file);
else
	ob.arg = file;
	wtfmex('gensys', ob.arg', sp_group_type);
	ob.nx = str2num(arg_get(ob.arg, 'nx'));
	ob.ny = str2num(arg_get(ob.arg, 'ny'));
	ob.nb = str2num(arg_get(ob.arg, 'nb'));
	ob.na = str2num(arg_get(ob.arg, 'na'));
end

ob.dims = [ob.nb * ob.na, ob.nx * ob.ny];

ob = class(ob, 'Gtomo2_wtfmex');

tmp = sum(ob) > 0;
ob.mask = reshape(tmp, ob.nx, ob.ny);

if ob.chat
	wtfmex('print'); % display header info
end
