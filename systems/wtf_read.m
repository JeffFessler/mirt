  function [mat nx ny nb na] = wtf_read(file, varargin)
%|function [mat nx ny nb na] = wtf_read(file, [options])
%| read aspire sparse matrix file (usually file.wtf)
%| and return matlab sparse matrix
%| option
%|	'parse'	0|1	if 1, do it by old matlab parsing way (old wtfread.m)
%| Copyright 2008-9-26, Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if nargin == 1 && streq(file, 'test'), wtf_read_test, return, end

arg.chat = 0;
arg.parse = 0;
arg = vararg_pair(arg, varargin);

if ~has_mex_jf
	fail('cannot read .wtf due to mex problem')
end

if arg.parse
	[mat nx ny nb na mask] = wtf_read_parse(file, arg.chat);
return
end

[buff nx ny nb na] = wtfmex('asp:read', file, int32(arg.chat));
[mat nx ny nb na is_tranpose] = wtfmex('asp:mat', buff, int32(arg.chat));
if is_tranpose
	mat = mat';
end


%
% wtf_read_parse()
% reads in a sparse matrix from Aspire format
% A is the returned sparse matrix
% n is a structure containing dimensions
% mask is the nx by ny support mask
%

function [A, nx, ny, nb, na, mask, n] = wtf_read_parse(file, chat)

machine = 'native';		% this should work
%machine = 'ieee-be';		% but try this if it doesn't...
fp = fopen(file, 'r', machine);
if (fp == -1), error fopen, end

%
% read until two form feeds
%
while (1)
	c = fread(fp, 1, 'char');
	if chat
		fprintf(1, '%c', char(c))	% echo header
	end
	if isempty(c)
		error eof
	end
	if (c == 12)	% form feed
		c == fread(fp, 1, 'char');
		if (c ~= 12)
			error only one form feed?
		else
			break;
		end
	end
end

%
% read binary header, extract some dimensions from it
%
tmp = fread(fp, 128/4, 'int');
type.group = tmp(1);
if (type.group ~= 1), error(sprintf('type.group = %d', type.group)), end
type.index = tmp(2);
if (type.index ~= 0), error 'type.index', end
type.value = tmp(3);
if (type.value ~= 0), error 'type.value', end
n.wt = tmp(6);
nx = tmp(12);
ny = tmp(13);
nb = tmp(14);
na = tmp(15);
if chat
	printf('nxy=%d,%d nba=%d,%d nwt=%d', ...
		nx, ny, nb, na, n.wt)
end

%
% read support
%
mask = fread(fp, nx*ny, 'uchar');
mask = reshape(mask, nx, ny);
if chat
	imagesc(mask');
end

%
% read length and offset
%
length = fread(fp, nx*ny, 'uint32');
offset = fread(fp, nx*ny, 'uint32');
if any(diff(offset) ~= length(1:end-1)), error bug, end

%
% read index and value arrays
%
index = fread(fp, n.wt, 'uint16');
value = fread(fp, n.wt, 'float32');

c = fread(fp, 1, 'char');
if ~isempty(c), error 'extra stuff in file?', end

%
% close file
%
if (fclose(fp)), error fclose, end

%
% generate sparse matrix
%
i = 1 + index;
j = zeros(n.wt,1);
n.col = nx * ny;
h = waitbar(0, 'Sparse matrix formation');
for ii=1:n.col
	t = [1:length(ii)] + offset(ii);
	j(t) = ii;
	waitbar(ii/n.col)
end
close(h)
if any(i <= 0) || any(j <= 0) || any(i > nb*na) || any(j > n.col)
	error 'bad indices'
end
if chat
	printf('wt value range [%g,%g]', min(value), max(value))
end
A = sparse(i, j, value, nb*na, n.col);


%
% wtf_read_test()
%
function wtf_read_test
nx = 6;
ny = 4;
nb = 8;
na = 5;
rng(0)
A = rand(nb*na, nx*ny);
A = double(single(A));

file = [test_dir 't.wtf'];
for row_grouped = 0:1
	pr row_grouped
	delete(file)
	wtf_write(file, A, nx, ny, nb, na, ...
		'chat', 0, 'row_grouped', row_grouped)
%	eval(['!wt test ' file])

	if 0 % test old "load" use
		[B mx my mb ma is_transpose] = wtfmex('load', file);
		pr is_transpose
		if is_transpose, B = B'; end
		jf_equal(A, B)
		jf_equal([mx my mb ma], [nx ny nb na])
	end

	if 1 % test new "asp:load" use
		[B mx my mb ma is_transpose] = wtfmex('asp:load', file);
		pr is_transpose
		if is_transpose, B = B'; end
		jf_equal(A, B)
		jf_equal([mx my mb ma], [nx ny nb na])
	end

	if 1 % test wtf_read
		[B mx my mb ma] = wtf_read(file);
		jf_equal(A, B)
		jf_equal([mx my mb ma], [nx ny nb na])
	end
end

ig = image_geom('nx', 22, 'ny', 20, 'dx', 2);
sg = sino_geom('par', 'nb', 24, 'na', 18, 'dr', 1.8);
ig.mask = ig.circ(ig.fov/2) > 0;

list = {'col', 'row'};
for ii=1:2
	gtype = list{ii};
	tmp = aspire_pair(sg, ig, 'support', 'array');
	buff = wtfmex('asp:gensys', tmp', gtype, uint8(ig.mask), int32(0));
	[mat mx my mb ma is_transpose] = wtfmex('asp:mat', buff, int32(0));
	jf_equal([mx my mb ma], [ig.nx ig.ny sg.nb sg.na])
	pr is_transpose
	if is_transpose, mat = mat'; end
	jf_equal(size(mat), [sg.nb*sg.na ig.nx*ig.ny])
end

% todo: test parse version?
