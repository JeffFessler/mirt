  function ob = vertcat(varargin)
%|function ob = vertcat(varargin)
%| called for ob = [a; b; ...] where any of them is a fatrix2

ob = cat(1, varargin{:});
