  function ob = horzcat(varargin)
%|function ob = horzcat(varargin)
%| called for ob = [a, b, ...] where any of them is a fatrix2

ob = cat(2, varargin{:});
