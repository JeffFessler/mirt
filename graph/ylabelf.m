 function h = ylabelf(varargin)
%function h = ylabelf(varargin)
%|
%| version of ylabel with built-in sprintf
%| and defaults to latex interpreter
%| also supports default font size from ir_fontsize()
%| Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

opt = {'fontsize', ir_fontsize('label')};
%opt = {opt{:}, 'fontname', 'times'};

varargin{1} = ir_strrep_tex(varargin{1}); % trick

if ir_is_octave
	tex = {};
	str = varargin{1};
	str = strrep(str, '$', ''); % remove $ needed for latex
	varargin{1} = str;
else
	tex = {'interpreter', 'latex'};
end

if 1
	str = varargin{1};
	str = strrep(str, '\', '\\'); % because of sprintf below
	% trick to replace "|" with "$|$" in strings
	if ~isempty(strfind(str, '|')) && isempty(strfind(str, '$'))
		str = strrep(str, '|', '$|$');
	end
	varargin{1} = str;
end

if im
	tmp = sprintf(varargin{:});
	tmp = strrep(tmp, '%', '\%'); % else latex will ignore %
	hh = ylabel(tmp, tex{:}, opt{:});
else
	hh = [];
end

if nargout
	h = hh;
end
