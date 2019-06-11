 function out = ir_fontsize(varargin)
%function out = ir_fontsize(varargin)
%|
%| set / return default font size for labels, title, etc
%|	'label'		return / set default for label
%|	'title'		return / set default for title
%|	'im_axes'	return / set default for im() axes
%|	'text'		return / set default for text
%|	'tick'		return / set default for xtick and ytick
%|
%|	'reset'		reset all values to defaults
%|
%| Copyright 2012-12-30, Jeff Fessler, University of Michigan


% ir_fontsize_reset
function fs = ir_fontsize_reset
	fs.label = 15;
	fs.title = 15;
	fs.text = 12;
	fs.im_axes = 11;
	fs.tick = 11;
end % ir_fontsize_reset


persistent fs % store preference
if ~isvar('fs') || isempty(fs)
	fs = ir_fontsize_reset;
end

if ~nargin, help(mfilename), pr fs, error(mfilename), end

switch nargin
case 1
	if streq(varargin{1}, 'reset')
		fs = ir_fontsize_reset;
		if nargout, out = fs; end
	else
		out = fs.(varargin{1});
	end

case 2
	if streq(varargin{1}, 'set') % set them all to same value
		val = varargin{2};
		if ischar(val)
			val = str2num(val);
		end
		fs.label = val;
		fs.title = val;
		fs.text = val;

	else
		val = varargin{2};
		if ischar(val)
			val = str2num(val);
		end

		fld = varargin{1};
		if ~isfield(fs, fld)
			fail('no field "%s"', fld)
		end
		fs.(fld) = val;
	end

otherwise
	fail('bug')
end

end % ir_fontsize
