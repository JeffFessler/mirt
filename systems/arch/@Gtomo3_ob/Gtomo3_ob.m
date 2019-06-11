 function ob = Gtomo3(sys_type, mask, nx, ny, nz, nthread, chat)
%function ob = Gtomo3(sys_type, mask, nx, ny, nz, nthread, chat)
% Construct Gtomo3 object, which can do 3d forward and backprojection
% using the ASPIRE cabilities for any system geometry.
% See Gtomo3_test.m for example usage.
% This object avoids the use of matlab's sparse matrices and instead
% uses the on-the-fly calculation of system matrix in f3d_mex.
% Basically, you create a f3d_mex system matrix by calling:
%		G = Gtomo3(sys_type, mask)
% and then you can use it thereafter by typing commands like
%		y = G * x;
% which will auto-magically call f3d_mex to do the multiplication
% in the compiled C mex program.
%
% in
%	sys_type	string	see 3D ASPIRE User's Guide
%	mask	[nx,ny,nz] or string
%		mask can be a filename, in which case the file is read.
%		mask can be empty (for some system types, including 3s),
%		in which case it is determined from backprojecting 1 vector.
%		Else, mask must be [nx,ny,nz] binary array of object support.
%	nthread			1, unless you have a two-processor system
%
% For more help, type 'f3d_mex'
%
% Copyright 01-04-22, Jeff Fessler, The University of Michigan

if nargin < 7, chat = 1; end
if nargin < 6, nthread = 1; end

%
% create default object, as required by Mathworks
%
ob = Gtomo2(chat);
ob.nz = 0;
ob = rmfield(ob, 'nb');
ob = rmfield(ob, 'na');
ob.n1 = 0;	% projection dimensions
ob.n2 = 0;
ob.n3 = 0;
ob.nthread = nthread;
ob.sys_type = '';
ob.file_mask = '';	% fix: delete

if nargin < 5	% required by Mathworks
	help(mfilename)
	warning 'Gtomo3 called with too few arguments!?'
	ob = class(ob, 'Gtomo3_ob');
	return
end

ob.nx = nx;
ob.ny = ny;
ob.nz = nz;
ob.sys_type = sys_type;

%
% if no mask, then backproject to find it!
%
if isempty(mask)
	f3d_mex('init', ob.sys_type, uint8(ones(nx,ny,nz)), ...
		int32(ob.nx), int32(ob.ny), int32(ob.nz), ...
		int32(ob.nthread), int32(ob.chat))
	warning('creating mask by summing will be a little slow')
	tmp = f3d_mex('dims')';
	mask = f3d_mex('back', ones(double(tmp), 'single'), int32(0)) > 0;
	f3d_mex('free')
end

if ischar(mask)
	mask = fld_read(mask); % read mask file
end

f3d_mex('init', ob.sys_type, uint8(mask), ...
	int32(ob.nx), int32(ob.ny), int32(ob.nz), ...
	int32(ob.nthread), int32(ob.chat))
tmp = double(f3d_mex('dims'));
[ob.n1, ob.n2, ob.n3] = deal(tmp(1), tmp(2), tmp(3));
ob.dims = [ob.n1 * ob.n2 * ob.n3, ob.nx * ob.ny * ob.nz];

%ob.is.empty = false;
ob = class(ob, 'Gtomo3_ob');

% check mask
ob.mask = mask;
if 1
	warning('checking mask, which may be slow and unnecessary')
	tmp = reshape(sum(ob) > 0, [ob.nx ob.ny ob.nz]);
	if ndims(tmp) ~= ndims(ob.mask), error 'wrong mask dims', end
	if any(size(tmp) ~= size(ob.mask)), error 'wrong mask size', end
	if any(tmp(:) ~= ob.mask(:))
		warning 'mask bug?  (could be due to poor sampling)'
		figure(1), clf, im(tmp, 'computed mask')
		figure(2), clf, im(ob.mask, 'input mask')
%		printf('range(tmp,2)'), disp(range(ob.mask,2)')
		disp('examine mask images then hit any key to continue:')
		pause
	end
	disp('mask checking completed')
end

if ob.chat
	f3d_mex('show', int32(ob.chat));	% display header info
end
