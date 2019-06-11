 function [wt scale180] = fbp_fan_short_wt(sg, varargin)
%function [wt scale180] = fbp_fan_short_wt(sg, [options])
%|
%| Sinogram weighting for fan-beam short scan, aka 'Parker weighting'
%|
%| in
%|	sg	strum		sino_geom (sg.orbit_start is ignored)
%|				or ct_geom
%|
%| option
%|	'type'	'parker'	from parker:82:oss (Med Phys 1982)
%|
%| out
%|	wt	[nb na]		parker weights; "excess" views have wt=0
%|	scale180		scale factor needed for back-projector	
%|
%| Usually you must multiply short fan-beam sinogram by
%| *both* the Parker weights "wt" and the scale factor
%| to get the proper results with my back projectors because
%| those back-projectors scale by the orbit and number of view angles.
%|
%| Copyright 2009-12-10, Jeff Fessler and Janghwan Cho, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(sg, 'test')
	[wt scale180] = fbp_fan_short_wt_test;
	if ~nargout, clear wt scale180, end
return
end

arg.type = 'parker';
arg = vararg_pair(arg, varargin);

switch arg.type
case 'parker'
	[wt scale180] = fbp_fan_short_wt_parker(sg);
otherwise
	fail('unknown type %s', arg.type)
end


% fbp_fan_short_wt_parker()
function [wt scale180] = fbp_fan_short_wt_parker(sg)

if sg.orbit < sg.orbit_short - eps
	warn('orbit %g is less than a short scan %g', sg.orbit, sg.orbit_short)
end

if sg.orbit > sg.orbit_short + eps
	warn('orbit %g exeeds short scan %g by %g views', ...
		sg.orbit, sg.orbit_short, ...
		(sg.orbit - sg.orbit_short) / (sg.ad(2) - sg.ad(1)));
end

nb = sg.nb;
na = sg.na;
bet = sg.ar - sg.ar(1); % trick: force 0 start, so this ignores orbit_start!
gam = sg.gamma;
[gg bb] = ndgrid(gam, bet);
gam_max = sg.gamma_max; % half of fan angle

fun = @(x) sin(pi/2 * x).^2; % smooth out [0,1] ramp
% todo: could use integrals of this function over the
% tiny angular range of each projection view so that
% the sum over beta of these functions is a constant.
% or use a bspline (quadratic?) that would be easier to integrate?

wt = zeros(nb,na); % any extra views will get 0 weight
ii = bb < 2 * (gam_max - gg); % 0 <= bb not needed
tmp = bb(ii) ./ (2 * (gam_max - gg(ii)));
wt(ii) = fun(tmp);

ii = 2 * (gam_max - gg) <= bb & bb < pi - 2 * gg;
wt(ii) = 1;

ii = pi - 2 * gg < bb & bb <= pi + 2 * gam_max;
tmp = (pi + 2*gam_max - bb(ii)) ./ (2 * (gam_max + gg(ii)));
wt(ii) = fun(tmp);

scale180 = sg.orbit / 180; % scale factor due to orbit and na in backprojector


% fbp_fan_short_wt_test
function [wt scale180] = fbp_fan_short_wt_test

% generate object, image geometry
down = 1;
ig = image_geom('nx', 512, 'dx', 1, 'down', down);
ell = [	0 0 200 160 0 1000;
	-70 0 50 50 0 1000;
	70 -0 2 2 0 1000];
xtrue = ellipse_im(ig, ell, 'oversample', 4);

orbit_start = 97; % stress test

%orbit_test = [0 20]; % for jang hwan cho
orbit_test = 0;

for io=1:numel(orbit_test)

	% sinogram geometry and fbp for full scan
	sg_360 = sino_geom('ge1', ...
		'orbit_start', orbit_start + orbit_test(io), 'down', down);
	sino_360 = ellipse_sino(sg_360, ell, 'oversample', 2);
	fbp_geom_360 = fbp2(sg_360, ig); % FBP reconstructed images
	fbp_w_360{io} = fbp2(sino_360, fbp_geom_360);

	sg_short = sino_geom('ge1', 'orbit', 'short', 'down', 4, ...
		'orbit_start', sg_360.orbit_start + orbit_test(io));
	if 1 % examine effect of extra views
		dd = sg_short.ad(2) - sg_short.ad(1);
		extra = 10;
		sg_short.na = sg_short.na + extra;
		sg_short.orbit = sg_short.orbit + extra * dd;
	end

	% sinogram of the object for each scan
	sino_short = ellipse_sino(sg_short, ell, 'oversample', 2);

	% apply parker weighting
	[wt scale180] = fbp_fan_short_wt(sg_short); % parker weighting
	sino_parker = sino_short .* wt;

	fbp_geom_short = fbp2(sg_short, ig); % FBP reconstructed images

	fbp_w_short = fbp2(sino_short, fbp_geom_short);
	fbp_w_parker{io} = scale180 * fbp2(sino_parker, fbp_geom_short);

	% plot
	im plc 4 3
	im subplot [1 2 3]
	im(rad2deg(sg_short.gamma), sg_short.ad, wt), cbar
	xlabel 'gamma [degrees]'
	ylabel 'beta [degrees]'
	title 'Parker weighting'
%	clim = [0 2000];
	clim = [800 1200];
	im(4, fbp_w_360{io}, clim, 'FBP: full scan'), cbar
	im(5, fbp_w_short, clim, 'FBP: short scan w/o parker weighting'), cbar
	im(6, fbp_w_parker{io}, clim, 'FBP, short scan w/ parker weighting'), cbar
	im(7, xtrue, clim, 'True'), cbar
	im subplot [8 9]
	plot([fbp_w_360{io}(:,end/2) fbp_w_short(:,end/2) fbp_w_parker{io}(:,end/2)])
	title 'middle slice, y = end/2'
	legend('full', 'short w/o parker', 'short w/ parker')

	im(10, sino_360)
	im(11, sino_short)
	im(12, sino_parker)
	if io < numel(orbit_test), prompt, end
end

if numel(orbit_test) > 1
	im_toggle(fbp_w_360{:}, [800 1200])
	prompt
	im_toggle(fbp_w_parker{:}, [800 1200])
	prompt
	im plc 1 2
	clim = [-1 1] * 100;
	im(1, fbp_w_360{2} - fbp_w_360{1}, clim), cbar
	im(2, fbp_w_parker{2} - fbp_w_parker{1}, clim), cbar
	keyboard
end

%clf, im_toggle(fbp_w_360 .* ig.circ, fbp_w_parker .* ig.circ)
%clf, plot(scale180 * mean(wt'), '.') % within 0.005 of 1

%yaxis_pi '0 p'
%plot(diff(wt,1))
% ir_savefig fig_tomo_fan_short_wt
