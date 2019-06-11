 function ob = Gtomo2_dsc(file, nthread, chat)
%function ob = Gtomo2_dsc(file, nthread, chat)
% Construct Gtomo2_dsc object, which can do 2d forward and backprojection
% using wtfmex's dsc,proj and dsc,back for any system geometry.
% See Gtomo2_dsc_test.m for example usage.
% This object computes the system matrix on the fly by calling the
% C-routines project_dsc and back_project_dsc
% Basically, you create a system matrix object by calling:
%		G = Gtomo2_dsc(file_dsc)
% and then you can use it thereafter by typing commands like
%		y = G * x;
% which will auto-magically call wtfmex to do the multiplication
% in the compiled C mex program.
% The ".dsc" file format is described in ASPIRE user's guide,
% a technical report available on my web page.
%
% Instead of a filename, one can instead pass a collection of
% argument pairs created using my arg_pair() routine.
%
% Use nthread = 2 for dual processor workstations
%
% Besides simple utilities like display, there are the following
% capabilities of this object:
%	y = G * x		forward projection
%	x = G' * y		back projection
%	y = Gt(:,ii)' * x	forward projection of a subset
%	x = Gt(:,ii) * y	back projection of a subset
%		the latter two only work if ii is a matrix of the form
%		ii = outer_sum(1:nb,(ia-1)*nb);
%
% Copyright 2001-08-30, Jeff Fessler, The University of Michigan

if ~isvar('chat') | isempty(chat)
	chat = 0;
end
if ~isvar('nthread') | isempty(nthread)
	nthread = 1;
end

%
% default object
%
ob = Gtomo2(chat);
ob.arg		= [];	% default dsc arguments
ob.nthread	= nthread;

if ~nargin | nargin > 3
	warning 'Gtomo2_dsc called with wrong number of arguments!?'
	help(mfilename)
	ob = class(ob, 'Gtomo2_dsc');
	return
end


%
% arg can be specified
%
if size(file,1) > 1
	ob.arg = file;

%
% or read from .dsc file
%
else
	fp = fopen(file, 'r');
	tmp = fgets(fp);
	while (tmp ~= -1)
		ob.arg = strvcat(ob.arg, tmp);
		tmp = fgets(fp);
	end
	fclose(fp);
end

ob.arg = remove_spaces(ob.arg);

ob.nx = str2num(arg_get(ob.arg, 'nx'));
ob.ny = str2num(arg_get(ob.arg, 'ny'));
ob.nb = str2num(arg_get(ob.arg, 'nb'));
ob.na = str2num(arg_get(ob.arg, 'na'));
ob.dims = [ob.nb*ob.na, ob.nx*ob.ny];

ob.mask = logical(wtfmex('dsc,mask', ob.arg', ob.chat));

ob = class(ob, 'Gtomo2_dsc');
