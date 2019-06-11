 function ir_im_pixel(x)
%function ir_im_pixel(x)
%|
%| Display image x in window with no menu or decoration and natural pixel size
%|
%| Jeff Fessler, based on:
%| http://stackoverflow.com/questions/6082896/remove-titlebar-from-matlab-gui-for-full-screen-display

if nargin < 1, help(mfilename), error(mfilename), end

switch class(x)
case 'uint8'
case {'double', 'single'}
	xmin = min(x(:));
	xmax = max(x(:));
	x = (x - xmin) / (xmax - xmin);
otherwise
	fail 'not done'
end

jimg = im2java(x);
frame = javax.swing.JFrame;
frame.setUndecorated(true)
icon = javax.swing.ImageIcon(jimg);
label = javax.swing.JLabel(icon);
frame.getContentPane.add(label);
frame.pack

nx_min = 100;

frame.setSize(max(nx_min, size(x,1)), max(nx_min, size(x,2)))

if 0
	screenSize = get(0,'ScreenSize');  %# Get the screen size from the root object
	frame.setSize(screenSize(3),screenSize(4));
	frame.setLocation(0,0);
end
frame.show
