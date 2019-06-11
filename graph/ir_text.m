 function h = ir_text(varargin)
%function h = ir_text(varargin)
%|
%| calls text() but specifies fontsize using ir_fontsize('text')
%|
%| 2013-09-08, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

if 1 % do TeX replacement in the first string argument
	for ii = 1:numel(varargin)
		if ischar(varargin{ii})
			varargin{ii} = ir_strrep_tex(varargin{ii}); % trick
			break
		end
	end
end

if 1
	tex = {'interpreter', 'latex'};
%	str = varargin{3};
%	str = strrep(str, '\', '\\'); % because of sprintf below
%	varargin{3} = str;
	varargin = {varargin{:}, tex{:}};
end

h = text(varargin{:}, 'fontsize', ir_fontsize('text'));
if ~nargout
	clear h
end
