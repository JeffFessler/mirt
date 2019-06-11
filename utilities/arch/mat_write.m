 function mat_write(file, data, varargin)
%function mat_write(file, data,
%		['-check', '-nocheck', '-type', datatype, '-name', varname])
%
% Save matlab "data" to "file", first checking if "file" exists!
% If file suffix is .fld, then save as AVS .fld file in xdr_float.
%
% datatype can be 'xdr_float' 'xdr_int' ...
% Jeff Fessler, The University of Michigan

if nargin == 1 && streq(file, 'test'), mat_write_test, return, end
if nargin < 2, ir_usage, end

% use .fld file if filename ends in .fld, or if filename is t0, t1, ...
isfld = 0;
if length(file) > 4
	isfld = streq(file((end-3):end), '.fld');
elseif length(file) == 2 && file(1) == 't'
	isfld = true;
end

check1st = true;
if isfld
	varname = '';
	datatype = 'xdr_float';
else
	varname = 'data';
	datatype = '';
end

while length(varargin)
	if streq(varargin{1}, '-nocheck')
		check1st = 0;
		varargin = {varargin{2:end}};
	elseif streq(varargin{1}, '-check')
		check1st = 1;
		varargin = {varargin{2:end}};
	elseif streq(varargin{1}, '-type')
		if length(varargin) < 2
			error '-type needs arg'
		else
			datatype = varargin{2};
			varargin = {varargin{3:end}};
		end
	elseif streq(varargin{1}, '-name')
		if length(varargin) < 2
			error '-name needs arg'
		else
			varname = varargin{2};
			varargin = {varargin{3:end}};
		end
	else
		error(sprintf('bad arg "%s"', varargin{1}))
	end
end

if isfld && ~isempty(varname)
	warning('variable name for .fld file ignored')
end
if ~isfld && ~isempty(datatype)
	warning('datatype for .mat file ignored')
end


if 2==exist(file) && check1st
	if has_aspire
		if isfld
			os_run(['op range ' file])
		else
			os_run(['op -chat 0 matls ' file])
		end
	end
	t = sprintf('file "%s" exists.  overwrite? [y|n]: ', file);
	t = input(t, 's');
	if ~streq(t, 'y')
		return
	end
end


% avs file
if isfld
	if check1st
		fld_write(file, data, '-check', datatype)
	else
		fld_write(file, data, '-nocheck', datatype)
	end


% matlab file
else	
	error 'writing to matlab file for reading by aspire no longer supported because mathworks changes the file format too often.  use fld_write'
	if ~streq(varname, 'data')
		eval([varname '= data'])	% caution: dangerous!
	end

	t = version;
	if t(1) == '4'
		save(file, varname)
	elseif t(1) == '5' || t(1) == '6'
		save(file, varname, '-v4')
	else
		error version
	end
end

printm('file "%s" written ok', file)


% mat_write_test()
function mat_write_test
dat = [5:8 1:4];
file = 'tmp.fld';
mat_write(file, dat, '-check', '-name', 'tester')
%d = load(file);
d = fld_read(file);
if any(d.tester ~= dat), error 'bug', else, printf('%s test ok', mfilename), end
delete(file)
