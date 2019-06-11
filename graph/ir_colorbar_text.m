 function ir_colorbar_text(str)
%function ir_colorbar_text(str)
%|
%| Add a colorbar and label it with string str.
%|
%| Needed as of octave 3.8.2 and matlab R2014b because they have incompatible
%| graphics support for colorbars.

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(str, 'test'), ir_colorbar_text_test, return, end

if ~im, return, end

h = colorbar;
if ir_is_octave
	subplot(h); ylabel(str)
else
	h.Label.String = str;
end


function ir_colorbar_text_test
im clf
im(eye(8))
ir_colorbar_text('test1')
