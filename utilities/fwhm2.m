  function [fw, angle, rad] = fwhm2(psf, varargin)
%|function [fw, angle, rad] = fwhm2(psf, [options])
%|
%| compute 2d fwhm of point-spread function using contourc()
%|
%| in
%|	psf [nx ny]	2D psf image
%|
%| option
%|	'dx'		pixel size (default 1)
%|	'dy'		pixel size (default dx)
%|	'level'		default 0.5 (for "half" maximum)
%|	'type'		'max' uses location of maximum (default)
%|			'centroid' takes 2D centroid around that max
%|			'user' uses user-defined center
%|	'center'	[cx cy] user-defined
%|	'chat'
%|
%| out
%|	fw		FWHM
%|
%| Copyright 2003-11-29, Jeff Fessler, University of Michigan
%| 2010-04-01 repaired bad bug in 'centroid' version that was former default
%| 2015-08-13 refinements by Valur Olafsson to better handle some PSFs

if nargin < 1, ir_usage, end
if streq(psf, 'test'), fwhm2_test, return, end

arg.dx = 1;
arg.dy = [];
arg.chat = (nargout == 0);
arg.level = 0.5;
arg.type = 'max';
arg.center = []; % for user-specified peak

%if nargin > 1 && isnumeric(varargin{1}) % old style: dx, chat
%	arg.dx = varargin{1}; varargin = {varargin{2:end}};
%	if length(varargin) > 1
%		arg.chat = varargin{2}; varargin = {varargin{2:end}};
%	end
%end

arg = vararg_pair(arg, varargin);
if isempty(arg.dy), arg.dy = arg.dx; end

if ~isreal(psf), warn 'unsure what happens for complex psf!', end

if min(size(psf)) < 11
	psf = padn(psf, max(size(psf), 23));
end
[nx ny] = size(psf);
ii = imax(psf, 2); % image maximum
cx = ii(1);
cy = ii(2);

switch arg.type
case 'max'
	% done

case 'user'
	if length(arg.center) ~= 2, error 'center should be [cx cy]', end
	cx = arg.center(1);
	cy = arg.center(2);

case 'centroid' % find center estimate by local centroid
	ix = [-5:5]';
	iy = [-5:5]';
	[iix iiy] = ndgrid(ix, iy);
	t = psf(cx + iix, cy + iiy);
	if arg.chat
		im(ix, iy, t)
	end
	ox = sum(t(:) .* iix(:)) / sum(t(:));
	oy = sum(t(:) .* iiy(:)) / sum(t(:));
	cx = cx + ox;
	cy = cy + oy;
	if ox ~= 0 || oy ~= 0
		printm('centroid offset %g %g', ox, oy)
	end

otherwise
	fail('bad type %s', arg.type)
end

if any(isnan(psf)), error 'psf contains nan', end
if min(psf(:)) > 1/2 * max(psf(:))
	warn 'psf does not decrease below max/2; could extrapolate?'
	fw = inf;
return
end

psf = double(psf); % stupid matlab
%cc = contourc(psf, [1e30 arg.level * max(psf(:))]);
cc = contourc(psf, [1e30 arg.level * psf(cx,cy)]); % 2015-08-13 change by valur olafsson
if isempty(cc), error 'empty contour?  check minimum!', end

% find the contour that has center of mass closest to center of PSF
if 1 % 2015-08-13 change by valur olafsson
	cc_colii = find(cc(1,:) == arg.level * psf(cx,cy));
	min_cc_mean = inf;
	min_cc_ii = 1;
	for ii = 1:length(cc_colii)
		cc_ii = cc_colii(ii);
		cc_mean = mean(cc(:,cc_ii+[1:cc(2,cc_ii)]), 2); % contour center of mass
		if sum(abs(cc_mean - [cy;cx])) < min_cc_mean
			min_cc_mean = sum(abs(cc_mean - [cy;cx]));
			min_cc_ii = cc_ii;
		end
	end
	cc = cc(:,min_cc_ii+[0:cc(2,min_cc_ii)]);
end

cc = cc(:,2:length(cc))';
cc = cc(:,[2 1]); % swap row,col or x,y

% check center pixel found
if arg.chat && im
	im plc 1 2
	im(1, psf)
	hold on
	plot(cc(:,1), cc(:,2), '+')
	plot(cx, cy, 'rx')
	titlef('length(cc)=%d', length(cc))
	hold off
end

xx = arg.dx * (cc(:,1) - cx); % physical coordinates
yy = arg.dy * (cc(:,2) - cy);
rsamp = sqrt(xx.^2 + yy.^2);
tsamp = atan2(yy,xx);
tsamp = rad2deg(tsamp);
[tsamp order] = sort(tsamp); % trick: monotone
rsamp = rsamp(order);
if abs(tsamp(end) - 180) < eps % -180 should match 180
	tsamp = [-180; tsamp];
	rsamp = [rsamp(end); rsamp];
end

angle = [0:180]'; % interpolate to equally spaced angles
% todo: replace with something better than linear interpolation
% todo: and make it handle 360 wrap around better?
r1 = interp1x(tsamp, rsamp, angle);
r2 = interp1x(tsamp, rsamp, angle-180);
rad = r1 + r2;
fw = mean(rad); % should the "average" be done s.t. for an ellipse contour
		% we get the avg of the two major axes?

if arg.chat && im
%	plot(tsamp, rsamp, 'o', angle, r1, '-', angle-180, r2, '-'), prompot
	im subplot 2
	plot(angle, rad, '-o', tsamp, 2*rsamp, 'yx')
	xlabel 'Angle \phi [degrees]'
	ylabel 'FWHM(\phi)'
	axisx(0,180), xtick([0 90 180])
	titlef('overall fwhm=%g', fw)
end


% fwhm2_test
function fwhm2_test
dx = 3;
dy = 2;
nx = 100;
ny = 80;
x = [-nx/2:nx/2-1]' * dx;
y = [-ny/2:ny/2-1]' * dy;
fx = 30;
fy = 16;
sx = fx / sqrt(log(256));
sy = fy / sqrt(log(256));
[xx yy] = ndgrid(x,y);
psf = exp(-((xx/sx).^2 + (yy/sy).^2)/2);
if im
	im pl 1 2
	im(1, x, y, psf, 'psf'), axis equal, axis image
end
[fw ang rad] = fwhm2(psf, 'dx', dx, 'dy', dy, 'chat', 1);
