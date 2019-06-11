 function x = fbp_dsc(sino, kernel, arg, varargin)
%function x = fbp_dsc(sino, kernel, arg, [options])
% Generic parallel-beam FBP reconstruction.
% Uses Gtomo2_dsc "system 9" backprojector (specific for pixel-driven FBP) 
% in:
%	sino	[nb,na] sinogram
%	kernel	for filtering
%	arg	argument array for Gtomo2_dsc, built using aspire_pair
%		or .dsc filename
% options
%	'nthread'
%	'chat'
% out:
%	x [np]	image estimate
%
% Copyright Apr 2000	Jeff Fessler	The University of Michigan

warning 'fbp_dsc is obsolete.  use fbp2 instead.  see fbp2_example'

if nargin == 1 && streq(sino, 'test'), fbp_dsc_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

opt.chat = 0;
opt.nthread = 1;
opt = vararg_pair(opt, varargin);

if ~isvar('kernel') || isempty(kernel)
	kernel = ones(3,1)/3;
end

sys = deblank(arg_get(arg, 'system', '0'));
if ~streq(sys, '9')
	printf('system type = "%s"', sys)
	warning 'only system 9 is likely to work'
end

G = Gtomo2_dscmex(arg, 'nthread', opt.nthread, 'chat', opt.chat);
mask = G.arg.mask;

if ~isempty(kernel) && any(size(kernel) ~= 1)
	sino = conv2(sino, kernel, 'same');	% filter (with zero pad)
end

sino = sino_ramp_filter(sino);
[nb na] = size(sino);
x = embed(G' * sino(:), mask);


%
% built-in test/example routine
% fbp_dsc_test()
%
function fbp_dsc_test
ig = image_geom('nx', 126, 'ny', 128, 'dx', 3);
ig.mask = ig.circ > 0;
sg = sino_geom('par', 'nb', 128, 'na', 140, 'dr', ig.dx, 'offset_r', 0.25);
arg2 = aspire_pair(sg, ig, 'system', 2, 'support', 'ellipse 0 0 50 50', ...
	'strip_width', sg.dr);
arg9 = aspire_pair(sg, ig, 'system', 9, 'support', 'ellipse 0 0 60 60');
G2 = Gtomo2_dscmex(arg2);
x = ig.circ(ig.fov/3) + ...
	ellipse_im(ig, [ig.nx/5 ig.ny/7 5 7 0 1], 'oversample', 3);
sino = G2 * x;
%sino = ellipse_sino(sg, );
fbp = fbp_dsc(sino, 1, arg9);
im pl 2 3
im(1, x, 'x'), cbar
im(4, sino, 'sinogram'), cbar
im(2, fbp, 'FBP (ramp)'), cbar
im(5, (fbp-x) .* ig.mask, 'error'), cbar
ix = 1:ig.nx; iy = ig.ny/2;
subplot(133), plot(ix, x(ix,iy), 'y-', ix, fbp(ix,iy), 'c-')
axis tight, legend('x', 'fbp')

if has_aspire
	dir = test_dir;
	f.sino = [dir 'sino.fld'];
	f.image = [dir 'image.fld'];
	f.dsc = [dir 't.dsc'];
	fld_write(f.sino, sino)
	char_array_write(f.dsc, arg2)
	f.win = 'boxcar,1,0,1';
	com = sprintf('echo y | i -chat 99 fbp dsc %s %s %s %s', ...
		f.image, f.sino, f.dsc, f.win);
	os_run(com)

	im_asp = fld_read(f.image);
	im(3, im_asp, 'aspire'), cbar
	good = x ~= 0;
	diff = im_asp - fbp;
	diff = good .* diff;
	diff(~good) = mean(diff(good));
	im(6, diff .* ig.mask, 'aspire-matlab'), cbar
%	im(4, (im_asp ~= 0) - (fbp ~= 0), 'support aspire-matlab'), cbar
end
