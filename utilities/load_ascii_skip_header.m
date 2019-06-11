  function data = load_ascii_skip_header(file, varargin)
%|function data = load_ascii_skip_header(file)
%| read ascii file, skipping lines that start with '#'
%| output is double array of numeric values
%| option
%|	'char'	0|1	if 1 then output cell array of each line
%| Copyright 2007, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
arg.char = false;
arg = vararg_pair(arg, varargin);

fid = fopen(file, 'r');
if fid == -1
	fail('problem opening "%s"', file)
end

if arg.char
	data = {};
else
	data = [];
end

while 1
	line = fgetl(fid);
	if ~ischar(line)
		if isempty(data)
			fail('file %s ended too soon', file)
		else
			break % go to fclose
		end
	end
	if ~length(line) || line(1) == '#'
		continue
	end
	if arg.char
		data{end+1} = line;
	else
		data(end+1,:) = sscanf(line, '%f');
	end
end

tmp = fclose(fid);
if tmp ~= 0
	fail('problem with fclose for "%s"', file)
end
