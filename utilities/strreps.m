function s = strreps(s, varargin)
%function s = strreps(s, f1, r1, f2, r2, ...)
% in string s, replace all f1 with r1, all f2 with r2, etc.

while(length(varargin))
	s = strrep(s, varargin{1}, varargin{2});
	varargin = {varargin{3:end}};
end
