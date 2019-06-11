 function A = Gtomo2_strip(sg, ig, varargin)
%function A = Gtomo2_strip(sg, ig, options)
%|
%| Generate a 2D system matrix for tomographic image reconstruction,
%| based on a square pixel basis and strip-integral detector model,
%| i.e., a rectangular detector response.
%| Works for parallel-beam, arc-fan, and flat-fan geometries.
%|
%| Closely matches Aspire "system 13" and "system 14".
%| The main limitation here is that it uses Matlab "sparse" data type,
%| so you will run out of memory for large image or sinogram sizes.
%|
%| in
%|	sg		sino_geom()
%|	ig		image_geom()
%|
%| options
%|	'strip_width'	if 0, then line integrals (default: sg.d)
%|	'single'	if 1, then double(single()) the values, cf aspire
%|	'gam_max'	max acceptance angle of detector [degrees] (def: 90)
%|
%| out:
%|	A [nb*na np]	Gsparse object, where np = sum(ig.mask(:))
%|
%| Copyright 2005-8-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sg, 'test'), Gtomo2_strip_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.gam_max = 90; % degrees
arg.strip_width = [];
arg.single = false;
arg.chat = false;
arg.class = 'fatrix2';
arg = vararg_pair(arg, varargin);

if isempty(arg.strip_width), arg.strip_width = sg.d; end

switch sg.type
case 'par'
	A = Gtomo2_strip_par(sg.nb, sg.na, sg.dr, sg.offset_r, ...
		sg.orbit, sg.orbit_start, arg.strip_width, ...
		ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
		ig.mask, arg.single);

case 'fan'
	A = Gtomo2_strip_fan(sg.nb, sg.na, sg.ds, sg.offset_s, ...
		sg.ar, arg.strip_width, ...
		deg2rad(arg.gam_max), ...
		sg.dsd, sg.dso, sg.dfs, ...
		sg.source_offset, ...
		ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
		ig.mask, arg.single, arg.chat);

otherwise
	fail('unknown geometry: %s', sg.type)
end

switch arg.class
case 'Fatrix'
	A = Gsparse(A, 'mask', ig.mask, 'idim', [ig.nx ig.ny], 'odim', [sg.nb sg.na]);
case 'fatrix2'
	A = Gmatrix(A(:,ig.mask(:)), 'imask', ig.mask, 'odim', [sg.nb sg.na]);
end


% Gtomo2_strip_par()
function A = Gtomo2_strip_par(nb, na, dr, offset_r, ...
	orbit, orbit_start, strip_width, ...
	nx, ny, dx, dy, offset_x, offset_y, mask, is_single);

% pixel centers
wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
wb = (nb-1)/2 + offset_r;

x = dx * ([0:nx-1] - wx);
y = dy * ([0:ny-1] - wy); % caution, may not match aspire if offset_y != 0
[x y] = ndgrid(x, y);
x = x(mask(:));
y = y(mask(:));
np = length(x);		% sum(mask(:)) - total # of support pixels

angle = deg2rad(orbit_start + [0:na-1]'/na * orbit);
cang = cos(angle);
sang = sin(angle);
tau = cang * x' + sang * y'; % [na np] projected pixel center
tau = tau / dr;
tau = tau + wb;

d_max = (abs(cang) + abs(sang)) / 2;
tau_max = (abs(dx) * d_max + strip_width/2) / dr;
ib_min = 1 + floor(tau - repmat(tau_max, 1, np));
M = ceil((abs(dx) * sqrt(2) + strip_width) / dr);
jj = find(mask(:))'; % all-column A
jj = repmat(jj, na, 1);
jj = col(jj); % so that na=1 case works
list.ii = [];
list.jj = [];
list.ss = [];
for mm=0:M-1
	ticker(mfilename, mm+1, M)
	ib = ib_min + mm;
	good = (ib >= 0) & (ib < nb);
	good = good(:);
	val = square_strip_int(dr * (ib - tau), ...
		repmat(angle, 1, np), 'dx', dx, 'sw', strip_width);
	val = val(:);
	ii = col(ib + repmat([0:na-1]'*nb, 1, np)); % sinogram index
	list.ii = [list.ii; 1+ii(good)];
	list.jj = [list.jj; jj(good)];
	list.ss = [list.ss; val(good)];
end
if is_single
	list.ss = double(single(list.ss)); % stupid matlab insists on double
end
A = sparse(list.ii, list.jj, list.ss, nb*na, nx*ny, length(list.ss));


% Gtomo2_strip_fan()
function A = Gtomo2_strip_fan(nb, na, ds, offset_s, ...
	beta, strip_width, ...
	gam_max, ... % maximum angle gamma accepted by detector [radians]
	dsd, dso, dfs, roff, ...
	nx, ny, dx, dy, offset_x, offset_y, mask, is_single, chat);

% pixel centers
wx = (nx-1)/2 + offset_x;
wy = (ny-1)/2 + offset_y;
wb = (nb-1)/2 + offset_s;

x = dx * ([0:nx-1] - wx);
y = dy * ([0:ny-1] - wy); % caution, may not match aspire if offset_y != 0
[x y] = ndgrid(x, y);
x = x(mask(:)); % [np,1]
y = y(mask(:));
np = length(x); % sum(mask(:)) - total # of support pixels
na = length(beta);

cbet = cos(beta);
sbet = sin(beta);

% s0 is "s" value for center of each pixel
gam0 = atan2(cbet * x' + sbet * y' - roff, ...
	dso + sbet * x' - cbet * y'); % [na np]
clear cbet sbet
if dfs == 0 % 3rd gen (arc)
	s0 = dsd * gam0;
elseif isinf(dfs)
	s0 = dsd * tan(gam0);
else
	error 'unsupported dfs'
end

ang0 = repmat(beta, [1 np]) + gam0;
mag0 = dso * cos(gam0) - roff * sin(gam0) ...
	+ repmat(x', [na 1]) .* sin(ang0) - repmat(y', [na 1]) .* cos(ang0);
gamgood = abs(gam0) < gam_max;
if chat, printm('gam0: %g %g', rad2deg(minmax(gam0(gamgood)))), end
clear gam0
if dfs == 0 % 3rd gen (arc)
	mag0 = mag0 / dsd;
elseif isinf(dfs)
	mag0 = mag0 / dsd ./ (1 + (s0 ./ dsd).^2);
else
	error 'dfs bug'
end
mag0min = min(mag0(gamgood));
if isempty(mag0min)
	A = sparse([], [], [], nb*na, nx*ny, 0);
	return
end

tau = s0 / ds + wb; % [na np], unitless

M = ceil((abs(dx)/mag0min * sqrt(2) + strip_width) / ds); % conservative!
if chat, printm('M=%d, mag0min=%g', M, mag0min), end
if M > 100
	printm('Warn: M=%d too large? probably x-ray source is too oblique', M)
	printm('type "return" and hope for the best, but it may be *slow*')
	printm('recommend using "gam_max" to impose incidence angle constraint')
	keyboard
end

ib_min = 1 + floor(tau - M/2);
jj = find(mask(:))';	% all-column A
jj = repmat(jj, na, 1); % [na np]
jj = col(jj); % so that na=1 case works
list.ii = [];
list.jj = [];
list.ss = [];
for mm=0:M-1
	ticker(mfilename, mm+1, M)
	ib = ib_min + mm;
	val = mag0 .* square_strip_int(ds * (ib - tau), ...
		ang0, 'dx', dx ./ mag0, ...
		'sw', strip_width);
	val = val(:); % for na=1 case
%	val = square_strip_int(ds * (ib - tau) .* mag0, ...
%		ang0, 'dx', dx, 'sw', strip_width .* mag0); % same!
	ii = col(ib + repmat([0:na-1]'*nb, 1, np)); % sinogram index
%	if chat, printm('%d of %d zeros', sum(val(:)==0), length(val(:))), end
	good = (ib >= 0) & (ib < nb);
	good = good & gamgood;
	good = good(:);
	good = good & (val > 0); % remove extra zeros due to conservative M
	list.ii = [list.ii; 1+ii(good)];
	list.jj = [list.jj; jj(good)];
	list.ss = [list.ss; val(good)];
end
if is_single
	list.ss = double(single(list.ss)); % stupid matlab insists on double
end
A = sparse(list.ii, list.jj, list.ss, nb*na, nx*ny, length(list.ss));


% Gtomo2_strip_test
function Gtomo2_strip_test

ig = image_geom('nx', 512, 'ny', 480', 'fov', 500, 'down', 8, ...
	'dy', '1.0*dx'); % stress-test with dy != -dx
ig.mask = ig.circ > 0;

for it=1:3
	if it == 1 % par
		sg = sino_geom('par', 'nb', 600, 'na', 480, 'dr', 1.2, ...
			'down', ig.down);
	elseif it == 2 % arc
		sg = sino_geom('fan', 'nb', 888, 'na', 984, 'ds', 1.0, ...
			'down', ig.down, ...
			'offset_s', 0*0.25, 'source_offset', 0*3.0, ...
			'dsd', 949, 'dod', 408, 'dfs', 0);
	elseif it == 3 % flat
		sg = sino_geom('fan', 'nb', 988, 'na', 984, 'ds', 1.0, ...
			'down', ig.down, ...
			'offset_s', 0*0.25, 'source_offset', 0*3.0, ...
			'dsd', 949, 'dod', 408, 'dfs', inf);
	end

	[x ell] = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	ya = ellipse_sino(sg, ell, 'oversample', 4);

	A = Gtomo2_strip(sg, ig, 'chat', 0);

	yd = A * x;
	sino = sg.zeros; sino(sg.nb/2, 10) = 1;
	b1 = A' * sino;
	bu = A' * sg.ones;
	max_percent_diff(min(bu(ig.mask)), max(bu(ig.mask)), 'backproject ones')

	if im
		im plc 2 2
		im(1, x, 'test image')
		im(2, ya, 'sinogram ya'), cbar
		im(4, yd, 'sinogram yd'), cbar
		im(3, yd-ya, 'yd-ya'), cbar
		xlabelf('nrms %5.2f%%', 100*nrms(yd(:), ya(:)))
		if 0
			im(1, ig.mask, 'support mask')
			im(2, b1, 'backproject 1 ray'), cbar
			im(3, bu, 'backproject ones'), cbar
		end
	end

	% verify consistency with Gtomo2_wtmex (aspire)
	if 1 && has_aspire %& streq(sg.type, 'par')
	prompt
		Aw = Gtomo2_wtmex(sg, ig, 'pairs', {'strip_width', sg.d});
		cpu etic
		for ii=1:1
			yw = Aw * x;
		end
		cpu etoc 'wtmex time'
		cpu etic
		for ii=1:1
			ys = A * x;
		end
		cpu etoc 'strip time'
		max_percent_diff(yw, ys, 'sino Gtomo2_wtmex vs Gtomo2_strip')
		if im
			im plc 2 3
			im(1, ya, 'ya'), cbar
			im(2, ys, 'ys'), cbar
			im(3, yw, 'yw'), cbar
			im(4, ys-yw, 'ys-yw'), cbar
			xlabelf('nrms %.2g%%', 100*nrms(ys(:), yw(:)))
			im(5, ys-ya, 'ys-ya'), cbar
			xlabelf('nrms %5.2f%%', 100*nrms(ys(:), ya(:)))
			im(6, yw-ya, 'yw-ya'), cbar
			xlabelf('nrms %5.2f%%', 100*nrms(yw(:), ya(:)))
		prompt
		end
	end
end
