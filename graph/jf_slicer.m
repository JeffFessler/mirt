 function jf_slicer(data, varargin)
%function jf_slicer(data, [options])
%|
%| slice 3d data interactively (along 3rd dimension).
%| Uses scroll wheel to sweep through all slices.
%|
%| todo: hsv!
%|
%| in
%|	data	[nx ny nz]
%|
%| options
%|	clim	[1 2]		clim arg to im()
%|	iz	[1]		initial slice (default: nz/2+1)
%|
%| Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if streq(data, 'test'), jf_slicer_test, return, end

arg.clim = [];
arg.iz = [];
arg = vararg_pair(arg, varargin);
if ~isreal(data)
	printm 'warning: taking abs of complex data'
	data = abs(data);
end

if isempty(arg.clim)
	arg.clim = double(minmax(data)');
end

nz = size(data, 3);

if isempty(arg.iz)
	iz = ceil((nz+1)/2);
else
	iz = arg.iz;
end

clamp = @(iz) max(min(iz, nz), 1);

stuff.data = data;
stuff.arg = arg;

%jf_slicer_call(iz, stuff)
%drawnow
jf_slicer_show

set(gcf, 'WindowScrollWheelFcn', @jf_slicer_scroll)
set(gcf, 'WindowKeyPressFcn', @jf_slicer_key)

%h = jf_add_slider('callback', @jf_slicer_call, 'data', {stuff}, ...
%	'min', 1, 'max', nz, 'sliderstep', [1 1]/(nz-1), 'value', iz);


function jf_slicer_scroll(src, event)
	iz = iz + event.VerticalScrollCount;
	iz = clamp(iz);
	jf_slicer_show
end % jf_slicer_scroll

function jf_slicer_key(src, event)
	key = get(gcf, 'CurrentCharacter');
	key = 0 + key;
	switch key
	case {28, 31, 45} % left arrow, down arrow, -
		iz = iz - 1;
	case {29, 30, 43} % right arrow, up arrow, +
		iz = iz + 1;
	otherwise
		if 49 <= key && key <= 57
			iz = key - 48; % 1-9
		end
	end
	iz = clamp(iz);
	pr iz
	jf_slicer_show
end % jf_slicer_key

function jf_slicer_show
	im(data(:,:,iz), arg.clim), cbar
	xlabelf('%d / %d', iz, nz)
end % jf_slicer_show

end % jf_slicer


function jf_slicer_call(iz, stuff)
persistent iz_save
if isempty(iz_save)
	iz_save = -1;
end
iz = round(iz);
if iz ~= iz_save
	iz_save = iz;
	arg = stuff.arg;
	im(stuff.data(:,:,iz), arg.clim), cbar
end

end % jf_slicer_call


function jf_slicer_test
data = reshape(1:7*8*9, 7, 8, 9);
im plc 1 2
im subplot 2
if im
	jf_slicer(data, 'clim', [0 500])
end
end % jf_slicer_test
