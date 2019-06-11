 function ir_fwrite(file, array, type, varargin)
%function ir_fwrite(file, array, type, varargin)
%|
%| Write array of "type" to file

if nargin < 3, ir_usage, end

arg.endian = 'native';
arg.check = 'fail';

if exist(file, 'file')
	switch arg.check
	case 'warn'
		warn('file "%s" exists, over-writing', file)
	case 'fail'
		fail('file "%s" exists', file)
	otherwise
		fail 'bug'
	end
end

fid = fopen(file, 'w', arg.endian);
if fid == -1, fail('cannot open %s', file), end

count = fwrite(fid, array, type);
if count ~= numel(array)
	fail('fwrite count problem')
end

status = fclose(fid);
if status ~= 0
	fail('fclose problem')
end
