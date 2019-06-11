 function fulldir = path_find_dir(part)
%function fulldir = path_find_dir(part)
% find full path to a directory given only the trailing part of it,
% searching the current matlab path
p = path;
dirend = [0 strfind(p, pathsep)];
for ii=2:length(dirend)
	fulldir = p([(dirend(ii-1)+1):(dirend(ii)-1)]);
	if streq(part, fulldir([max(end-length(part)+1,1):end]))
		return
	end
end
error(sprintf('could not find directory fragment "%s" in path', part))
