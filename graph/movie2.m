 function [mov avi] = movie2(x, varargin)
%function [mov avi] = movie2(x, [options])
%|
%| in
%|	x	[nx ny nz]	sequence of 2d frames
%| option
%|	clim	"color" limits			default: minmax(x)
%|	file	name of avi file		default: [test_dir 'tmp.avi']
%|	cmap	colormap			default: gray(256)
%|	fps	frames per second in avifile	default: 15
%| out
%|	mov	movie object
%|	avi	avifile object
%|
%| make matlab movie from 3d array (e.g., iterations)
%| make avi movie from 3d array (e.g., iterations)
%|
%| If no output arguments, then display movie
%| Copyright Aug 2000, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), movie2_test, return, end

arg.clim = [];
arg.file = [test_dir 'tmp.avi'];
arg.cmap = gray(256);
arg.fps = 15; % default for avifile
arg = vararg_pair(arg, varargin);

if ~isreal(x)
	printm 'warning: abs of complex data'
	x = abs(x);
end

if isempty(arg.clim)
	arg.clim = minmax(x)';
end
pr arg.clim

%clf
%fig = figure;
%fig = gcf;
%set(fig, 'DoubleBuffer', 'on');
%set(gca, 'NextPlot', 'replace', 'Visible', 'off')
x = max(x, arg.clim(1));
x = min(x, arg.clim(2));
scale = 255. / single(arg.clim(2) - arg.clim(1));
x = scale * single(x - arg.clim(1)); % need single for uint16 inputs

% avi = avifile(arg.file, 'fps', arg.fps);
avi = VideoWriter(arg.file, 'Grayscale AVI');
% 'VideoCompressionMethod', 'None');
open(avi)
F.colormap = [];
nz = size(x,3);
clf
for iz=1:nz
	ticker(mfilename, iz, nz)
	t = x(:,:,iz)';
	t = uint8(t);

	mov(iz) = im2frame(t, arg.cmap); % for matlab "movie"

	if 0 % old way for avifile
		F.cdata = repmat(t, [1 1 3]);
		avi = addframe(avi,F);
	end

	im(t, arg.clim, sprintf('%d', iz-1)), drawnow

	if 0 && iz == 1
		set(gca, 'nextplot', 'replacechildren')
%		set(gcf, 'Renderer', 'zbuffer')
	end
%	frame = getframe;
	frame = mov(iz);
	frame.colormap = [];
%	pr iz
%	try
	writeVideo(avi, frame)
%	catch
%	printm 'error: writevideo'
%	keyboard
%	end

%	F = getframe(gca);
end
% avi = close(avi);
close(avi)

if ~nargout
	if im
		movie(mov, 2);
	end
	clear mov avi
end


%
% movie2_test()
%
function movie2_test
down = 30;
cg = ct_geom('fan', 'ns', round(888/down), 'nt', 64, 'na', round(984/down), ...
	'ds', 1.0*down, 'down', 1, ... % only downsample s and beta
	'offset_s', 0.25, ... % quarter detector
	'offset_t', 0.0, ...
	'ztrans', 1*300, ... % stress test with helix
	'dsd', 949, 'dod', 408, 'dfs', inf); % flat detector
% 'dsd', 949, 'dod', 408, 'dfs', 0); % 3rd gen CT

ell = [0*50 0*50 0*50 200 100 100 90 0 10];
proj = ellipsoid_proj(cg, ell, 'oversample', 2);

[mov avi] = movie2(proj);
if im
	movie(mov)
%	keyboard
end
