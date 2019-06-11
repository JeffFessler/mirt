  function wtf_write(file, mat, nx, ny, nb, na, varargin)
%|function wtf_write(file, mat, nx, ny, nb, na, [options])
%| write (usually sparse) matrix mat to file (usually file.wtf)
%|	mat can be Gsparse Fatrix or Gmatrix fatrix2 also.
%| option
%|	'row_grouped'	1|0	1 to write in row grouping (default: 0 col)
%|	'chat'
%! Copyright 2008-9-26, Jeff Fessler, University of Michigan

if nargin == 1 && streq(file, 'test'), wtf_write_test, return, end
if nargin < 6, help(mfilename), error(mfilename), end

arg.chat = 0;
arg.row_grouped = 0;
arg = vararg_pair(arg, varargin);

if ~has_mex_jf
	fail('cannot write .wtf due to mex problem')
end

switch class(mat)
case 'Fatrix'
	if ~isvar('mat.arg.G')
		fail 'unknown Fatrix'
	end
	mat = mat.arg.G; % for Gsparse object
	if ncol(mat) ~= nx * ny
		if ~isvar('mat.arg.mask')
			fail 'no mask?'
		end
		mat = embed(mat, mat.arg.mask); % todo: can it work?
	end

case 'fatrix2'
	if ~isvar('mat.arg.matrix')
		fail 'unknown Fatrix'
	end
	ob = mat;
	mat = ob.arg.matrix; % for Gmatrix object
	if ncol(mat) ~= nx * ny
		tmp = fatrix2_embed(ob.imask, ob.idim, mat');
		tmp = shiftdim(tmp, numel(ob.idim));
		tmp = reshape(tmp, nb*na, nx*ny);
		mat = tmp;
	end
end

if ~issparse(mat)
	mat = sparse(mat); % ensure sparse
end

jf_equal([nb*na nx*ny], size(mat)) % verify size

wtfmex('asp:save', file, mat, int32(nx), int32(ny), ...
	int32(nb), int32(na), int32(arg.chat), int32(arg.row_grouped));


% wtf_write_test
function wtf_write_test
run_mfile_local('wtf_read test')
