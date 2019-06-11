 function im_toggle(i1, i2, varargin)
%function im_toggle(i1, i2, [..., options for im()])
%|
%| toggle between two or more images via keypress (octave) or click (matlab)
%|
%| Jeff Fessler

if nargin == 1 && streq(i1, 'test'), im_toggle_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% find leading additional arguments corresponding to images
iall = {i1, i2};
names = {inputname(1), inputname(2)};
ii = 3;
while length(varargin) && isequal(size(i2), size(varargin{1}))
	iall = {iall{:}, varargin{1}};
	varargin = {varargin{2:end}};
	names{ii} = inputname(ii);
	ii = ii + 1;
end

for ii=1:length(names)
	names{ii} = sprintf(['toggle i%d: ' names{ii}], ii);
end

if length(varargin) && ( ...
	streq(varargin{1}, 'mip') || ...
	streq(varargin{1}, 'sum') || ...
	streq(varargin{1}, 'mid') )
	for ii=1:length(iall)
		iall{ii} = jf_mip3(iall{ii}, 'type', varargin{1});
	end
	varargin = {varargin{2:end}};
end

if ~im, return, end

% toggle between two or more images

ft.args = {{}, {'horizontalalignment', 'right'}};
ft.pos = [0.01, 0.99];

if ir_is_octave % && exist('uicontrol') ~= 2
	im_toggle_key(iall, ft, names, varargin);
else % matlab
	im_toggle_gui(iall, ft, names, varargin);
end


% im_toggle_gui()
function im_toggle_gui(iall, ft, names, args)

im clf
ht = uicontrol('style', 'text', 'string', 'test', ...
                'units', 'normalized', 'position', [0.1 0.03 0.8 0.03]);
data = {iall, ft, args, names, ht};
im_toggle_call(0, data{:}, true)

jf_add_slider('callback', @im_toggle_call, 'data', data, ...
	'sliderstep', [0.999 1] / (length(iall)-1));

%set(gcf, 'WindowScrollWheelFcn', @im_toggle_scroll) % todo
%set(gcf, 'WindowKeyPressFcn', @im_toggle_key) % todo


% im_toggle_key()
% old way based on key presses
function im_toggle_key(iall, ft, names, args)

while (1)
	for ii=1:length(iall)
%		im_toggle_call(ii/length(iall), iall, ft, args, names)
		im clf
		im(iall{ii}, args{:})
		fig_text(ft.pos(2-mod(ii,2)), 0.01, ...
			names{ii}, ...
			ft.args{2-mod(ii,2)})

%		pause
		in = input('hit enter for next image, or "q" to quit ', 's');
		if streq(in, 'q')
			set(gca, 'nextplot', 'replace')
			return
		end
	end
end


function im_toggle_call(value, iall, ft, im_args, names, ht, reset)
persistent counter
persistent last_value
if isempty(last_value)
	last_value = value;
end
if isempty(counter) || isvar('reset')
	counter = 0;
end
if value == last_value % click
	counter = 1 + counter;
	if counter > length(iall)
		counter = 1;
	end
else % slide
	counter = 1 + floor(0.9999 * length(iall) * value);
	last_value = value;
end

im(iall{counter}, im_args{:})
set(ht, 'string', names{counter})


function im_toggle_test
if 0 % 2d
	nx = 20;
	i1 = eye(20);
	i2 = flipud(i1);
	i3 = 1 - i1;
	im_toggle(i1, i2, i3, [0 2])
end

if 1 % 3d mip3
	ig = image_geom('nx', 2^5, 'ny', 2^5-1', 'nz', 2^5-3, 'dx', 1);
	t1 = ellipsoid_im(ig, [4 3 2 ig.fovs ./ [3 3 4] 0 0 1]);
	t2 = flipdim(t1, 1);
	t3 = flipdim(t1, 2);
	im_toggle(t1, t2, t3, 'sum')
end
