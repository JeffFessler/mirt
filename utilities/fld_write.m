 function fld_write(file, data, varargin)
%function fld_write(file, data, [options])
%|
%| write data into AVS format .fld file
%| see ASPIRE user's guide (or AVS online documentation) for file format.
%|
%| todo: generalize to include 'extent' options and default
%|
%| options
%|	'check'		0|1	1: report error if file exists (default: 0)
%|	'dir'		char	directory name to prepend file name
%|	'type'		char	e.g., 'short_be'.  default is 'xdr_float'
%|	'head'		strings	comment information for file header
%|	'raw'		0|1	1: put raw data in name.raw, header in name.fld
%|					where file = name.fld
%|	'endian'	str	'ieee-le' or 'ieee-be' for little/big endian
%|				default inferred from 'type'
%|				with bias to 'ieee-le' (changed 2016-06-20)
%|
%| Copyright 2003-5-13, Jeff Fessler, University of Michigan

%| options (obsolete way)
%|	'-check'		report error if file already exists
%|	'-nocheck'		overwrite file even if it already exists
%|	'header' [strings]	comment information for file header
%|	datatype		e.g., 'short_be'.  default is 'xdr_float'

if nargin == 1 && streq(file, 'test'), fld_write_test, return, end
if nargin < 2, ir_usage, end

if length(varargin)
	arg1 = varargin{1};
	if streq(arg1, '-check') || streq(arg1, '-nocheck') ...
		|| streq(arg1, 'header') || streq(arg1, 'byte')
		warn 'old style syntax; please upgrade!'
		fld_write_old(file, data, varargin{:});
	elseif streq(arg1, 'new')
		fail 'the ''new'' option no longer needed - remove!'
	else
		fld_write_new(file, data, varargin{:});
	end
else
	fld_write_new(file, data);
end


% fld_write_new()
function fld_write_new(file, data, varargin)
arg.check = 0;
arg.raw = 0;
arg.type = '';
arg.head = '';
arg.dir = '';
arg.endian = ''; % infer below from arg.type

arg = vararg_pair(arg, varargin);

if isempty(arg.endian)
	switch arg.type
	case {'short_be', 'int_be', 'float_be', 'double_be', ...
		'xdr_int', 'xdr_double' 'xdr_float'}
		arg.endian = 'ieee-be';
	otherwise
		arg.endian = 'ieee-le';
	end
end

if isempty(arg.type)
	switch arg.endian
	case 'ieee-le'
		arg.type = 'float_le';
	case 'ieee-be'
		arg.type = 'xdr_float';
	otherwise
		fail('bad endian %s', arg.endian)
	end
end

% try to check that endian and data type match
if streq(arg.type, 'xdr_', 4) && streq(arg.endian, 'ieee-le')
	warn('for xdr_ data types like "%s", endian should be ieee-be; switching', arg.type)
	arg.endian = 'ieee-be';
end
if any(strfind(arg.type, '_le')) && streq(arg.endian, 'ieee-be')
	warn('for _le data types, endian should be ieee-le')
	arg.endian = 'ieee-le';
end
if any(strfind(arg.type, '_be')) && streq(arg.endian, 'ieee-le')
	warn('for _be data types, endian should be ieee-be')
	arg.endian = 'ieee-be';
end

if ~isempty(arg.dir)
	file = [arg.dir filesep file];
	printm('file = "%s"', file)
end

fld_write_do(file, data, arg.check, arg.type, arg.endian, arg.head, arg.raw);


% fld_write_old()
function fld_write_old(file, data, varargin);

% defaults
check1st = 0;
datatype = 'xdr_float'; % default
endian = 'ieee-be';
header = '';
while length(varargin)
	arg1 = varargin{1};
	varargin = {varargin{2:end}};
	if streq(arg1, '-check')
		check1st = 1;
	elseif streq(arg1, '-nocheck')
		check1st = 0;
	elseif streq(arg1, 'header')
		header = varargin{1};
		varargin = {varargin{2:end}};
	else
		datatype = arg1;
	end
end

fld_write_do(file, data, check1st, datatype, endian, header, 0);


% fld_write_do()
function fld_write_do(file, data, check1st, datatype, endian, header, raw)

% check if file already exists
loop = check1st & exist(file, 'file');
while (loop)
	t = sprintf('file "%s" exists.  overwrite? [y|n|d]: ', file);
	t = input(t, 's');
	if streq(t, 'd')
		t = fld_read(file);
		printf('max %% diff = %g', max_percent_diff(t, data))
		continue
	end
	if streq(t, 'y')
		break
	else
		return
	end
end

if raw
	fileraw = file; fileraw(end+[-2:0]) = 'raw';
	loop = check1st && exist(fileraw);
	while (loop)
		t = sprintf('file "%s" exists.  overwrite? [y|(n)]: ', fileraw);
		t = input(t, 's');
		if streq(t, 'y')
			break
		else
			return
		end
	end
end

% data type
if 1
	switch datatype
	case {'byte', 'uint8'}
		datatype = 'byte';
		format = 'uint8';
	case {'short_be', 'short_sun'}
		datatype = 'short_be';
		format = 'int16';
	case {'short_le', 'short', 'int16'}
		datatype = 'short_le';
		format = 'int16';
	case 'float'
		format = 'float';
		endian = 'native'; % not portable
	case {'float_le', 'float'}
		format = 'float';
	case {'float_be', 'xdr_float'}
		datatype = 'xdr_float';
		format = 'float32';
	otherwise
		error(['datatype ' datatype ' not implemented yet'])
	end
end

% open avs file for writing
fid = fopen(file, 'w', endian);
if fid == -1, fail('cannot open %s', file), end
if raw
	fraw = fopen(fileraw, 'w', endian);
else
	fraw = fid;
end

% write header
ndim = ndims(data);

fprintf(fid, '# created by %s\n', mfilename);
if ~isempty(header)
	for ii=1:nrow(header)
		fprintf(fid, '# %s\n', header(ii,:));
	end
end
fprintf(fid, 'ndim=%d\n', ndim);
%fprintf(fid, 'nspace=%d\n', ndim);
for ii=1:ndim
	fprintf(fid, 'dim%d=%d\n', ii, size(data,ii));
end
fprintf(fid, 'data=%s\n', datatype);
fprintf(fid, 'veclen=1\n');
fprintf(fid, 'field=uniform\n');
if raw
	fprintf(fid, 'variable 1 file=%s filetype=binary\n', fileraw);
else
	fprintf(fid, '\f\f');	% two form feeds: char(12)
end

%
% finally, write the binary data
%

count = fwrite(fraw, data, format);
if count ~= numel(data)
	disp([count size(data)])
	fail('file count=%d vs data=%d', count, prod(dims))
end

if fclose(fid), error 'fclose error?', end
if raw
	if fclose(fraw), error 'fclose error?', end
end



function fld_write_test1(file, data, varargin)
fld_write(file, data, varargin{:})
tmp = fld_read(file);
if any(tmp ~= data), keyboard, error 'test failed', end
delete(file)


function fld_write_test
file = [test_dir 'tmp.fld'];
data = [5:8; 1:4];
%fld_write(file, data)
%fld_write(file, data, 'new', 'check', 1, 'endian', 'ieee-le')
if 1 % test raw/header
	fld_write_test1(file, data, 'check', 1, 'endian', 'ieee-le', 'raw', 1)
	delete([test_dir 'tmp.raw'])
end

formats = {'float_be', 'float_le', 'float', 'xdr_float'};
for ii=1:length(formats)
	format = formats{ii};
	pr format
	fld_write_test1(file, data, 'check', 1, 'type', format)
	%		 'endian', 'ieee-le', ... % nah, defer from type
end
	format = 'short_le';
	pr format
	fld_write_test1(file, int16(data), 'check', 1, 'type', format)

printm 'passed'
