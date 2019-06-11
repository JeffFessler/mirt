  function proj = square_strip_int(rr, angles, varargin)
%|function proj = square_strip_int(rr, angles, varargin)
%|
%| strip-integrals of a square pixel basis object model,
%| normalized such that the units are length, not area.
%| In other words, the detector response is h(r) = (1/w) rect(r/w).
%| in
%|	rr	[]	radial sample locations, relative to 0.
%|	angles	[]	radians
%|
%| options
%|	dx		pixel_size (default: 1)
%|	sw		strip_width (default: dx)
%|
%| Copyright 2005-8-1, Jeff Fessler, University of Michigan

if nargin == 1, square_strip_int_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.dx = 1; % pixel_size
arg.sw = []; % strip_width
arg = vararg_pair(arg, varargin);
if isempty(arg.sw), arg.sw = arg.dx; end

proj = square_strip_int_do(rr, angles, arg.dx, arg.sw);


%
% square_strip_int_do()
%
function proj = square_strip_int_do(rr, angles, dx, sw)

cang = abs(cos(angles));
sang = abs(sin(angles));
dmax = (cang + sang) / 2;
dbrk = abs(cang - sang) / 2;
lang = 1 ./ max(cang,sang);

% see book: s,geom,par2,strip
if sw == 0
	proj = zeros(size(rr));
	tt = rr ./ dx;
	ii = (-dmax < tt) & (tt < -dbrk);
	proj(ii) = (tt(ii) + pick(dmax, ii)) ...
		./ (pick(dmax, ii) - pick(dbrk, ii));
	ii = abs(tt) <= dbrk;
	proj(ii) = 1;
	ii = (dbrk < tt) & (tt < dmax);
	proj(ii) = (pick(dmax, ii) - tt(ii)) ...
		./ (pick(dmax, ii) - pick(dbrk, ii));

	proj = (dx .* lang) .* proj;

else
	proj = gam((rr-sw/2) ./ dx, (rr+sw/2) ./ dx, lang, dmax, dbrk);
	proj = dx.^2 ./ sw .* proj; % normalize by strip width, so length units
end


function out = pick(vv, ii)
if length(vv) > 1
	out = vv(ii);
else
	out = vv;
end

function g = gam(a, b, lang, dmax, dbrk)
g ...
= gam1(max(a, -dmax), min(b, -dbrk), lang, dmax, dbrk) ...
+ gam2(max(a, -dbrk), min(b, +dbrk), lang) ...
+ gam3(max(a, +dbrk), min(b, +dmax), lang, dmax, dbrk);

function g = gam1(a, b, lang, dmax, dbrk)
dmax(dmax == dbrk) = dmax(dmax == dbrk) + 10*eps; % trick
g = lang/2 ./ (dmax - dbrk) .* ((b + dmax).^2 - (a + dmax).^2) .* (b > a);

function g = gam2(a, b, lang)
g = lang .* (b - a) .* (b > a);

function g = gam3(a, b, lang, dmax, dbrk)
dmax(dmax == dbrk) = dmax(dmax == dbrk) + 10*eps; % trick
g = lang/2 ./ (dmax - dbrk) .* ((a - dmax).^2 - (b - dmax).^2) .* (b > a);


function square_strip_int_test

r = linspace(-3,3,4001);
a = [0 pi/8 pi/4 pi/3 pi/2];
[rr aa] = ndgrid(r, a);
dx = 2.0; sw = 0.5;
proj = square_strip_int(rr, aa, 'dx', dx, 'sw', sw);

if im
	clf, plot(r, proj, '-')
	title(sprintf('dx=%g sw=%g', dx, sw))
	axis([minmax(r)' 0 2.9])
	legend(cellstr(num2str(rad2deg(a'))))
end
