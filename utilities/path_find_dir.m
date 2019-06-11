 function fulldir = path_find_dir(part)
%function fulldir = path_find_dir(part)
%|
%| find full path to a directory given only the trailing part of it,
%| searching the current matlab path
%|
if nargin < 1, ir_usage, end

p = path; % long string separated by ':' (pathsep)
dirend = [0 strfind(p, pathsep)]; % all the ':' locations
for ii=2:numel(dirend)
	fulldir = p([(dirend(ii-1)+1):(dirend(ii)-1)]); % one dir at a time
	if streq(part, fulldir([max(end-numel(part)+1,1):end])) % match end
		return
	end
end
fail('could not find directory fragment "%s" in path', part)
