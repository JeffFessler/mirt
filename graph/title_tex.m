  function h = title_tex(arg, varargin)
%|function title_tex(arg, varargin)
%| show title with tex interpreter
%| also supports default font size from ir_fontsize()

if nargin < 1, help(mfilename), error(mfilename), end

opt = {'fontsize', ir_fontsize('title')};

if im
	h = title(arg, 'interpreter', 'tex', varargin{:}, opt{:});
else
	h = [];
end

if ~nargout
	clear h
end
