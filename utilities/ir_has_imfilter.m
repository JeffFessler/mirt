 function out = ir_has_imfilter 
%function out = ir_has_imfilter 
%| determine if we have the "real" imfilter from image processing toolbox or not

out = (exist('imfilter') == 2) && ir_has_image_processing_toolbox;


function out = ir_has_image_processing_toolbox 
v = ver;
if ~isstruct(v)
	fail 'cannot determine if this is octave or matlab!?'
end

% note: must use strcmp not my streq to support string-cell comparisons
tmp = strcmp('Image Processing Toolbox', {v.Name});
out = any(tmp);
