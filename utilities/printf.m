 function printf(str, varargin)
%function printf(str, varargin)
% shorthand for disp(sprintf())
% note: fprintf() is overloaded to work without a file id, so i could have
% used it instead.  oh well.
% Jeff Fessler
if nargin < 1, ir_usage, end
disp(sprintf(str, varargin{:}))
