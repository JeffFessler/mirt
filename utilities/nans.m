  function out = nans(varargin)
%|function out = nans(varargin)
%| single precision nan for array initialization
out = nan(varargin{:}, 'single');
