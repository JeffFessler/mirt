 function opt = vararg_pairs(opt, varargin)
%function opt = vararg_pairs(opt, varargin)
% obsolete function: use vararg_pair
if nargin < 1, ir_usage, end
warning 'vararg_pairs(opt, varargin{:}) is obsolete, use vararg_pair(opt, varargin) instead'

opt = vararg_pair(opt, varargin)
