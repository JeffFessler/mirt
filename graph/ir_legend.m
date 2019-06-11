 function h = ir_legend(varargin)
%function h = ir_legend(varargin)
%|
%| calls legend() but specifies fontsize using ir_fontsize('text')
%| and asks for latex interpreter
%|
%| 2013-09-08, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

if ~iscell(varargin{1})
	fail 'need cell argument'
end
if 1 % do TeX replacement in the cell (legend) arguments
	for ii = 1:numel(varargin{1})
		varargin{1}{ii} = ir_strrep_tex(varargin{1}{ii}); % trick
	end
end

if 1
	tex = {'interpreter', 'latex'};
	varargin = {varargin{:}, tex{:}};
end

%if ~contains({varargin{2:end}}, 'location') % 2016b
if ~any(strcmp({varargin{2:end}}, 'location'))
	varargin = {varargin{:}, 'location', 'best'};
end

h = legend(varargin{:}, 'fontsize', ir_fontsize('text'));
if ~nargout
	clear h
end
