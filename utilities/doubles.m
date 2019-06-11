 function d = doubles(s)
%function d = doubles(s)
%
% convert a value to a double, but check first that it is numeric!
%
% Copyright 2003-3-11	Jeff Fessler	The University of Michigan

if nargin < 1, ir_usage, end

if islogical(s)
	d = double(s);
elseif ~isnumeric(s)	% because matlab 6.5 thinks logical is not numeric!
	error 'not numeric'
else
	d = double(s);
end
