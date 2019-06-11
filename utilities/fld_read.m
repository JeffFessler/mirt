 function [data, coord, dims] = fld_read(file, varargin)
%function [data, coord, dims] = fld_read(file, options)
%| Read data from AVS format .fld file.
%|
%| options
%|	'dir'	''	prepend file name with this directory (default: '')
%|	'raw'	0|1	1: return raw data class (default), 0: return doubles
%|	'slice' int	specify which slice to read from 3D file. (0 ... nz-1)
%|	'chat'	0|1	enable verbosity
%|	'dim_only' 0|1	returns dims. data and coord are equal to [].
%|	'coord' 0|1	returns coordinates too (default: 0)
%|	'coord_format'	default: 'n'; see fopen() for machine format options
%|			(needed for some UM .fld files with 'vaxd' coordinates)
%| out
%|	data		stores the element values.
%|	coord		vector that contains the coordinate values in the order
%|				X axis, Y axis, Z axis, ...
%|	dims		vector that specifies the length of each dimension.
%|
%| see ASPIRE user's guide for documentation
%|
%| Copyright 2003-5-13, Jeff Fessler & Rongping Zeng, University of Michigan

if nargin == 1 && streq(file, 'test'), fld_write('test'), return, end
if nargin < 1, ir_usage, end

% defaults
arg.slice = [];
arg.dir = [];
arg.raw = true;
arg.chat = false;
arg.dim_only = false;
arg.coord = false;
arg.coord_format = 'n'; % default 'native
arg = vararg_pair(arg, varargin);

if ~isempty(arg.dir)
	file = [arg.dir filesep file];
end

fid = fopen(file, 'r');
if fid == -1, fail('cannot open %s', file), end
%onCleanup(@() fclose(fid)); % todo: use this?

formfeed = char(12);	% form feed
newline = sprintf('\n');
is_external_file = false;

% read header until we find the end of file or the 1st form feed
header = '';
while (1)
	[inchar count] = fread(fid, 1, '*char'); % read one character

	% end of file means external file (or error)
	if count ~= 1
		if ~isempty(strfind(header, 'file='))
			is_external_file = true;
			break		% end of header file!
		else
			disp(header)
			error 'end of file before form feeds?'
		end
	end

	% form feed means embedded file
	if inchar == formfeed
		[inchar count] = fread(fid, 1, '*char');
		if count ~= 1, error 'end of file before 2nd form feed?', end
		if inchar ~= formfeed, error 'not two form feeds?', end
		if (arg.chat)
			printm('embedded data file')
		end
		break
	end

	% otherwise append this character to header string
	header = [header inchar];
end, clear inchar formfeed count newline

header = string_to_array(header); % convert to array
if (arg.chat)
	disp(header)
end

% parse header to determine data dimensions and type
ndim = arg_get(header, 'ndim');
dims = zeros(1,ndim);
for ii=1:ndim
	dims(ii) = arg_get(header, sprintf('dim%d', ii));
end
fieldtype = arg_get(header, 'field', '%s');
datatype = arg_get(header, 'data', '%s');
if arg_get(header, 'veclen') ~= 1, error 'only veclen=1 done', end
if arg.chat
	printm('ndim=%d', ndim)
	printm('dim%d=%d ', [1:ndim; dims])
end

% process dim_only option
if arg.dim_only
	data = [];
	coord = [];
	fclose(fid);
return
end

% external file (binary data in another file)
% fix: external ASCII files to be implemented
skip = 0;
if is_external_file
	fclose(fid);
	extfile = arg_get(header, 'file', '%s');

	filetype = arg_get(header, 'filetype', '%s');
	if arg.chat, printm('Current file = "%s", External file = "%s", type="%s"', ...
			file, extfile, filetype), end

	if ~isempty(strfind(col(header')', 'skip='))
		skip = arg_get(col(header')', 'skip');
	end

	if ~streq(filetype, 'multi')
		if ~exist(extfile, 'file')
			fdir = file;
			slash = strfind(fdir, '/');
			if isempty(slash)
				fail('cannot find external file %s', extfile)
			end
			fdir = fdir(1:slash(end));
			extfile = [fdir extfile]; % add directory
			if ~exist(extfile, 'file')
				fail('no external ref file %s', extfile)
			end
		end
	else
		num_of_files = str2num(extfile);
		fdir = file;
		slash = strfind(fdir, '/');
		if isempty(slash)
			fail('cannot find external file %s', extfile)
		end
		fdir = fdir(1:slash(end));
	end
else
	filetype = '';
	extfile = '';
end

% finally, read the binary data
[format endian bytes] = datatype_fld_to_mat(datatype);

if streq(filetype, 'multi') % multi file reading
	if arg.coord, fail 'coord not implemented for multi files', end
	coord = [];
	data = fld_read_multi(file, fid, arg, ...
		dims, datatype, fieldtype, fdir, header, ...
		is_external_file, extfile, format, endian, bytes, skip);

else % single file reading
	[data, coord] = fld_read_single(file, fid, arg, ...
		dims, datatype, fieldtype, ...
		is_external_file, extfile, format, endian, bytes, skip);

	fclose(fid);
end


%
% fld_read_single()
%
function [data, coord] = fld_read_single(file, fid, arg, ...
		dims, datatype, fieldtype, ...
		is_external_file, extfile, format, endian, bytes, skip)

% reopen file to same position, with appropriate endian too.
if is_external_file
	if isempty(endian)
		fid = fopen(extfile, 'r');
	else
		fid = fopen(extfile, 'r', endian);
	end
	if fid == -1, fail('cannot open external %s', extfile), end
else
	position = ftell(fid);
	if fclose(fid) ~= 0, warning fclose, end
	if isempty(endian)
		fid = fopen(file, 'r');
	else
		fid = fopen(file, 'r', endian);
	end;
	if fid == -1, fail('cannot re-open %s', file), end
	if fseek(fid, position, 'bof') ~= 0, error fseek, end
end

if skip
	if fseek(fid, skip, 'cof') ~= 0, error fseek, end
end

% if a single slice to be read, then skip to it
if ~isempty(arg.slice)
	if length(dims) ~= 3
		if arg.slice == 0 || arg.slice == -1
			dims(3) = 1;
			arg.slice = 0;
		else
			error 'slice only good for 3d files', end
		end
	if (arg.slice < 0), arg.slice = arg.slice + dims(3); end % trick for negative slice
	if (arg.slice < 0 || arg.slice >= dims(3)), error 'bad slice', end
	offset = bytes * arg.slice * dims(1) * dims(2);
	if fseek(fid, offset, 'cof') ~= 0, error fseek, end
	rdims = dims(1:2);
else
	rdims = dims;
end

% read binary data (converting to double) and reshape appropriately
% to avoid converting to double, we would need to preface format with a '*',
% per fread documentation.

if arg.raw
	format = ['*' format];
end

[data count] = fread(fid, prod(rdims), format);
if count ~= prod(rdims)
	pr rdims
	fid
	fail('file count=%d vs data=%d', count, prod(rdims))
end
if length(rdims) > 1
	data = reshape(data, rdims);
end


% reopen file to the same position, but to specified machine format
% to read in coordinates.
% UM RadOnc .fld files sometimes use 'vaxd' format for coordinates,
% which are no longer supported by matlab, so they may not work.
if arg.coord
	position = ftell(fid);
	fclose(fid);
	if streq(datatype, 'xdr_short')
		fid = fopen(file, 'r', 'ieee-be');
	else
		fread_coord = @(file,n,type) fread(file,n,type);
		try
			fid = fopen(file, 'r', arg.coord_format);
		catch
			if ~streq(arg.coord_format, 'vaxd')
				fail('fopen(%s) with format "%s" failed', ...
					file, arg.coord_format)
			end
			if exist('freadVAXD') ~= 2
				fail(['need to get "freadVAXD from ', ...
	' http://www.mathworks.com/matlabcentral/fileexchange/22675'])
			end
			fid = fopen(file, 'r', 'ieee-le'); % trick
			fread_coord = @(f,n,type) freadVAXD(f,n,type); % trick
		end
	end
	fseek(fid, position, 'bof');
	switch fieldtype
	case 'uniform'
		% printm('Number of coordinate values is 0')
		coord = [];

	case 'rectilinear'
		cor_num = sum(dims);
		% trick: so single-slice works:
		fseek(fid, -(cor_num*4+1), 'eof');
		[coord count] = fread_coord(fid, inf, 'float');
		if count ~= cor_num
			error('error in reading coordinate values')
		end

	otherwise
		fail('fieldtype "%s" not implemented yet!', fieldtype)
	end
else
	coord = [];
end


%
% fld_read_multi()
% multiple file reading (someshs@umich.edu)
% reusing code from single external file/same file reading above
%
% first implemented for the case when separate slices are saved
% as different binary files.
%
% todo: implement skip
% todo: implement vaxd
%
function data = fld_read_multi(file, fid, arg, ...
		dims, datatype, fieldtype, fdir, header, ...
		is_external_file, extfile, format, endian, bytes, skip)

if ~is_external_file
	error 'file=num_of_files usage missing when filetype=multi is used'
end
if skip
	error 'skip not yet implemented for filetype=multi.'
end

extfile_list = '';

% if a line of fld file does not contain a '='
% then we use that line as one of the filenames
for ii=1:size(header,1)
	temp = header(ii,:);
	if isempty(strfind(temp,'='))
		extfile_list = strvcat(extfile_list,temp);
	end
end
if arg.chat
	printm('current dir = %s', fdir)
	printm 'file list:-'
	printm '=============================================================='
	extfile_list
	printm '=============================================================='
end

% setup format per fread documentation
if arg.raw
	format = ['*' format];
end

% if a single slice to be read, then read that particular slice file
if ~isempty(arg.slice)

	if length(dims) ~= 3
		if arg.slice == 0 || arg.slice == -1
			dims(3) = 1;
			arg.slice = 0;
		else
			error 'slice only good for 3d files', end
		end
	if (arg.slice < 0), arg.slice = arg.slice + dims(3); end % trick for negative slice
	if (arg.slice < 0 || arg.slice >= dims(3)), error 'bad slice', end

	extfile_slice = deblank([fdir extfile_list(arg.slice,:)]);
	fid = fopen(extfile_slice, 'r', endian);
	if fid == -1, fail('cannot open external %s', extfile_slice), end
	rdims = dims(1:2);

	% read in and reshape the data
	[data count] = fread(fid, prod(rdims), format);
	if count ~= prod(rdims)
		disp(rdims)
		fid
		fail('file count=%d vs data=%d', count, prod(rdims))
	end
	fclose(fid);

else
	rdims = dims;
	rdims_slice = dims(1:2);
	data = [];
	count = 0;
	ticker reset
	nfile = size(extfile_list,1);
	for ii=1:nfile % read in each file
		ticker(mfilename, ii, nfile)
		extfile_slice = deblank([fdir extfile_list(ii,:)]);
		fid = fopen(extfile_slice, 'r', endian);
		if fid == -1, fail('cannot open external %s', extfile_slice), end

		[data_new count_new] = fread(fid, prod(rdims_slice), format);
		data = [data; data_new];
		count = count + count_new;
		fclose(fid);
	end

	if count ~= prod(rdims)
		fail('file count=%d vs data=%d', count, prod(rdims))
	end
	if length(rdims) > 1
		data = reshape(data, rdims);
	end

end


%
% string_to_array()
% convert long string with embedded newlines into string array
%
function header = string_to_array(header_lines)
newline = sprintf('\n');

% ensure there is a newline at end, since dumb editors can forget...
if header_lines(end) ~= newline
	header_lines = [header_lines newline];
end

ends = strfind(header_lines, newline);
if length(ends) <= 0
	error 'no newlines?'
end

header = header_lines(1:(ends(1)-1));
for ll = 2:length(ends)
	line = header_lines((ends(ll-1)+1):(ends(ll)-1));
	header = strvcat(header, line);
end

% strip comments (lines that begin with #)
header(header(:,1) == '#',:) = [];


%
% arg_get()
% parse an argument from header, of the name=value form
%
function arg = arg_get(head, name, type)
if ~isvar('type')
	type = '%d';
end
for ll = 1:nrow(head)
	line = head(ll,:);
	start = strfind(line, [name '=']);
	if ~isempty(start)
		if length(start) > 1, error 'bug: multiples?', end
		line = line((start+length(name)+1):end);
		[arg, count, err] = sscanf(line, type, 1);
		if ~isempty(err)
			error(err)
		end
		return
	end
end
fail('could not find %s in header', name)


%
% datatype_fld_to_mat()
% determine matlab format from .fld header datatype
%
function [format, endian, bytes] = datatype_fld_to_mat(datatype)
switch datatype
case 'byte'
	format = 'uint8';
	endian = 'ieee-be'; % irrelevant
	bytes = 1;

case 'short'
	format = 'short';
	endian = ''; % native short - not portable
	bytes = 2;
case {'short_be', 'short_sun', 'xdr_short'}
	format = 'int16';
	endian = 'ieee-be';
	bytes = 2;
case 'short_le'
	format = 'int16';
	endian = 'ieee-le';
	bytes = 2;

case 'int'
	format = 'int';
	endian = ''; % native int - not portable
	bytes = 4;
case 'int_le'
	format = 'int32';
	endian = 'ieee-le';
	bytes = 4;
case {'int_be', 'xdr_int'}
	format = 'int32';
	endian = 'ieee-be';
	bytes = 4;

case 'float'
	format = 'float';
	endian = ''; % native float - not portable
	bytes = 4;
case 'float_le'
	format = 'float32';
	endian = 'ieee-le';
	bytes = 4;
case  {'float_be', 'xdr_float'}
	format = 'float32';
	endian = 'ieee-be';
	bytes = 4;

case 'double'
	format = 'double';
	endian = ''; % native double - not portable
	bytes = 8;
case 'double_le'
	format = 'double';
	endian = 'ieee-le';
	bytes = 8;
case {'double_be', 'xdr_double'}
	format = 'float64';
	endian = 'ieee-be';
	bytes = 8;
otherwise
	error(['format "' datatype '" not yet implemented. ask jeff!'])
end
