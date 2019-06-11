 function varargout = num2list(x)
%function varargout = num2list(x)
%
% convert a numeric vector input into a comma separated list of output values
% DOES NOT WORK!

if nargin < 1, num2list_test, ir_usage, end
%nargout
%nargout = numel(x)
%varargout = cell(size(x));

for ii=1:numel(x)
	varargout{ii} = x(ii);
end

function num2list_test
zeros(num2list([2 3 4 5]))
