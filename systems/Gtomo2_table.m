 function ob = Gtomo2_table(sg, ig, table, varargin)
%function ob = Gtomo2_table(sg, ig, table, options)
%|
%| Generate a 2D parallel-beam tomographic system model
%| based on finely tabulated footprints.
%|
%| in
%|	sg	strum		sino_geom()
%|	ig	strum		image_geom()
%|	table	[Ktab Ltab]	footprint table, see book: s,geom,par2
%|				or see below for pre-defined options
%| options
%|	'nthread'		# threads for parellizing; default: jf('ncore')
%|
%| out
%|	ob	[nb*na np]	Fatrix system "matrix", np = sum(ig.mask(:))
%|
%| pre-defined "options" for table:
%|	{'linear', 'Ltab', [1000]}
%|		pixel-driven linear interp, not recommended except for testing
%|	{'square/strip', 'Ltab', [1000], 'strip_width', [dr]}
%|		square-pixel, strip-integral model, supercedes "system 2"
%|	{'dd2', 'Ltab', [1000], 'strip_width', [dr]}
%|		approximation that is numerically equivalent to distance driven
%|		(not recommended except for exploring the properties of DD )
%|	{'la45', ...} local angle approximation every 45 degrees
%|
%|	the following make an exact table based on mojette radial sampling:
%|	{'mojette,back1'}	linear interpolation for backprojection
%|	{'mojette,linear'}
%|	{'mojette,square/strip'}
%|
%| Copyright 2005-7-31, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sg, 'test'), Gtomo2_table_test, return, end
if nargin < 3, ir_usage, end

if streq(sg.type, 'fan'), fail 'only parallel-beam and mojette done', end

% options
arg.sg = sg;
arg.ig = ig;
arg.table = table;
arg.chat = 0;
arg.class = 'fatrix2';
arg.nthread = jf('ncore');
arg = vararg_pair(arg, varargin);

arg = Gtomo2_table_setup(arg);

arg.nd = sg.nb * sg.na;
arg.np = sum(ig.mask(:));
arg.odim = [sg.nb sg.na];
dim = [arg.nd arg.np]; % trick: make it masked by default!

% build object
switch arg.class
case 'Fatrix'
	ob = Fatrix(dim, arg, 'caller', 'Gtomo2_table', ...
		'abs', @(ob) ob, ... % nonnegative
		'forw', @Gtomo2_table_forw, 'back', @Gtomo2_table_back, ...
		'mtimes_block', @Gtomo2_table_mtimes_block);

case 'fatrix2'
	forw = @(arg, x) Gtomo2_table_block_forw_array(arg, x, 1, 1);
	back = @(arg, y) Gtomo2_table_block_back_array(arg, y, 1, 1);
	forw_block = @(arg, x, iblock, nblock) ...
		Gtomo2_table_block_forw_array(arg, x, iblock, nblock);
	back_block = @(arg, y, iblock, nblock) ...
		Gtomo2_table_block_back_array(arg, y, iblock, nblock);
	ob = fatrix2('mask', ig.mask, 'arg', arg, 'does_many', 1, ...
		'abs', @(ob) ob, ... % nonnegative
		'forw', forw, 'back', back, ...
		'forw_block', forw_block, 'back_block', back_block);

otherwise
	fail('bad class "%s"', arg.class)
end


% Gtomo2_table_setup()
function arg = Gtomo2_table_setup(arg)

sg = arg.sg;
ig = arg.ig;
arg.mask8 = uint8(ig.mask);
arg.angles = sg.ar;

arg.ptab_arg = {int32(arg.nthread), int32(sg.nb), sg.d, ...
	sg.offset, ig.dx, ig.dy, ig.offset_x, ig.offset_y};

if ~has_mex_jf, fail 'jf_mex required', end

if isnumeric(arg.table)
	arg.table_type = 'numeric';
elseif iscell(arg.table)
	if streq(arg.table{1}, 'mojette', 7)
		arg = Gtomo2_table_setup_mojette(arg);
		arg.proj_str = 'moj2,proj';
		arg.back_str = 'moj2,back';
		arg.ptab_arg = {arg.ptab_arg{[1:2 4:8]}}; % no dr
	else
		arg = Gtomo2_table_setup_cell(arg);
		arg.proj_str = 'ptab,par2,proj';
		arg.back_str = 'ptab,par2,back';
	end
else
	fail 'unknown table class'
end

arg.table = single(arg.table);


% Gtomo2_table_setup_cell()
function arg = Gtomo2_table_setup_cell(arg)

arg.table_type = arg.table{1};
arg.table = {arg.table{2:end}};
opt.chat = 0;
opt.Ltab = 1000;

switch arg.table_type

case 'linear' % pixel-driven linear interpolation
	opt.Ktab = 2;
	opt = vararg_pair(opt, arg.table);
	tri = @(t) (1-abs(t)) .* (abs(t) < 1);
	ir_min = @(t) 1 + floor(t - 1);
	[kk ll] = ndgrid(0:(opt.Ktab-1), 0:opt.Ltab-1);
	tau = (ll + 1/2) / opt.Ltab;
	arg.table = arg.ig.dx * tri(ir_min(tau) + kk - tau);

	if im && opt.chat
		plot(ll(1,:), arg.table', '.'), grid
	end

case {'square/strip', 'dd2', 'la45'}
	opt.Ktab = [];
	opt.strip_width = [];
	opt = vararg_pair(opt, arg.table);
	if isempty(opt.strip_width)
		opt.strip_width = arg.sg.d;
	end
	tmax = (arg.ig.dx/sqrt(2) + opt.strip_width/2) / arg.sg.dr;
	if isempty(opt.Ktab), opt.Ktab = ceil(2*tmax); end
	ir_min = @(t, tau_max) 1 + floor(t - tau_max);
	[kk ll ang] = ndgrid(0:opt.Ktab-1, 0:opt.Ltab-1, arg.angles);
	tau = (ll + 1/2) / opt.Ltab;

	t = arg.sg.dr * (ir_min(tau, opt.Ktab/2) + kk - tau);
	switch arg.table_type
	case 'dd2' % distance driven equivalent
		mang = max(abs(cos(ang)), abs(sin(ang)));
		arg.table = square_strip_int(t, 0*ang, ...
			'dx', arg.ig.dx * mang, ...
			'sw', opt.strip_width) ./ mang.^2;

	case 'la45' % local angle approximation, every 45 degrees
		ang8 = floor((ang + pi/8) / (pi/4)) * pi/4; % nearest 45
%		pr unique(rad2deg(ang8))
%		clf, plot(rad2deg(ang(:)), rad2deg(ang8(:))), prompt
		mang = max(abs(cos(ang)), abs(sin(ang))) ...
			./ max(abs(cos(ang8)), abs(sin(ang8)));

		if 1
			ia = (ang8 == pi/4);
			mang(ia) = 1; % todo: why better?
			%mang(ia) = 1 ./ mang(ia);
			%mang = 1;
			%mang = 1 ./ mang;
		end

		arg.table = square_strip_int(t, ang8, ...
			'dx', arg.ig.dx * mang, ...
			'sw', opt.strip_width) ./ mang.^2;

	case 'square/strip'
		arg.table = square_strip_int(t, ang, ...
			'dx', arg.ig.dx, ...
			'sw', opt.strip_width);
	otherwise
		fail 'bug'
	end

	if im && opt.chat
		plot(tau(1,:,1), arg.table(:,:,1)', '.'), grid
	end

otherwise
	fail('unknown table type %s', arg.table_type)
end

% function for displaying the table...
arg.t_fun4 = @(k, l, K, L) 1 + floor((l+1/2)/L - K/2) + k - (l+1/2)/L;
arg.t_fun = @(k, l) arg.t_fun4(k, l, opt.Ktab, opt.Ltab);
arg.tab_opt = opt;


% Gtomo2_table_setup_mojette()
function arg = Gtomo2_table_setup_mojette(arg)

arg.table_type = arg.table{1};
arg.table = {arg.table{2:end}};
sg = arg.sg;
ig = arg.ig;

if sg.dx ~= abs(ig.dx) || abs(ig.dx) ~= abs(ig.dy)
	fail('mojette dx=%g vs dx,dy=%g,%g', sg.dx, ig.dx, ig.dy)
end

% pixel-driven linear interpolation (intended for backprojection only)
switch arg.table_type
case 'mojette,back1'
	Ktab = 2;
	[rtab, ang, dr] = mojette_rad_ang(ig.nx, ig.ny, ig.dx, ig.dy, ...
		ig.offset_x, ig.offset_y, ...
		sg.nb, sg.offset_r, arg.angles, Ktab);
	tri = @(t) (1-abs(t)) .* (abs(t) < 1);
	arg.table = tri(rtab ./ dr); % [0,1] linear interpolation coefficients
%	clf, plot(reshape(rtab./dr, [], arg.na), reshape(arg.table, [], arg.na), '.'), prompt

% forward projector that relates to pixel-driven linear interpolation
case 'mojette,linear'
	Ktab = 2;

	[rtab ang dr] = mojette_rad_ang(ig.nx, ig.ny, ig.dx, ig.dy, ...
		ig.offset_x, ig.offset_y, ...
		sg.nb, sg.offset_r, arg.angles, Ktab);

	tri = @(t) (1-abs(t)) .* (abs(t) < 1);
	arg.table = ig.dx^2 ./ dr .* tri(rtab ./ dr);


% square pixel / strip integral model, strip_width = dr = dx max(|cos|,|sin|)
case 'mojette,square/strip'
	opt.strip_width = [];
	opt = vararg_pair(opt, arg.table);
	if isempty(opt.strip_width) || opt.strip_width > 0
		Ktab = 3; % special case due to strip width choice
	else
		Ktab = 2; % special case due to strip width choice
	end

	[rtab, ang] = mojette_rad_ang(ig.nx, ig.ny, ig.dx, ig.dy, ...
		ig.offset_x, ig.offset_y, ...
		sg.nb, sg.offset_r, arg.angles, Ktab);

	if Ktab == 3
		sw = abs(ig.dx) * max(abs(cos(ang)), abs(sin(ang)));
	else
		sw = 0;
	end

	arg.table = square_strip_int(rtab, ang, 'dx', ig.dx, 'sw', sw);

%	clf, plot(reshape(rtab, [], arg.na), reshape(arg.table, [], arg.na), '.')

otherwise
	fail 'unknown mojette table type'
end


% mojette_rad_ang()
% radial and angular samples for complete mojette table
% out
%	rtab	[Ktab max(nx,ny) na]
%	ang	""
%	dr	""
%	irmin	""
%
function [rtab ang dr irmin] = mojette_rad_ang(nx, ny, ...
	dx, dy, offset_x, offset_y, ...
	nb, offset_r, angles, Ktab)

[kk mn ang] = ndgrid(0:Ktab-1, 0:max(nx,ny)-1, angles);
% ir_min = @(t, tau_max) 1 + floor(t - tau_max);
wr = (nb - 1)/2 + offset_r;
wx = (nx - 1)/2 + offset_x;
wy = (ny - 1)/2 + offset_y;

cang = cos(ang);
sang = sin(ang);
cang(abs(ang - pi/2) < eps) = 0;
sang(abs(ang - pi/2) < eps) = 1;

tau = zeros(size(kk));

ia = abs(cang) >= abs(sang);

ii = ia & (dx * cang > 0);
tau(ii) = wr - wx + (mn(ii) - wy) .* ...
	(dy * sang(ii)) ./ (dx * cang(ii));

ii = ia & (dx * cang < 0);
tau(ii) = wr - (nx-1-wx) + (mn(ii) - wy) .* ...
	(dy .* sang(ii)) ./ abs(dx * cang(ii));

ia = abs(cang) < abs(sang);

ii = ia & (dy * sang > 0);
tau(ii) = wr - wy + (mn(ii) - wx) .* ...
	(dx * cang(ii)) ./ (dy * sang(ii));

ii = ia & (dy * sang < 0);
tau(ii) = wr - (ny-1-wy) + (mn(ii) - wx) .* ...
	(dx * cang(ii)) ./ abs(dy * sang(ii));

irmin = 1 + floor(tau - Ktab/2);

if any(irmin(:) < 0) || any(irmin(:) > nb-Ktab)
	minmax(irmin)
	fail 'truncated mojette sinogram.  increase nb'
end
dr = max(abs(dx * cos(ang)), abs(dy * sin(ang)));
rtab = (irmin + kk - tau) .* dr;


% Gtomo2_table_forw(): y = G * x
function y = Gtomo2_table_forw(arg, x)
y = Gtomo2_table_mtimes_block(arg, 0, x, 1, 1);


% Gtomo2_table_back(): x = G' * y
function x = Gtomo2_table_back(arg, y)
x = Gtomo2_table_mtimes_block(arg, 1, y, 1, 1);


% Gtomo2_table_mtimes_block()
function y = Gtomo2_table_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	y = Gtomo2_table_block_back(arg, x, istart, nblock);
else
	y = Gtomo2_table_block_forw(arg, x, istart, nblock);
end


% Gtomo2_table_block_forw()
function y = Gtomo2_table_block_forw(arg, x, istart, nblock)

[x ei] = embed_in(x, arg.ig.mask, arg.ig.np);
y = Gtomo2_table_block_forw_array(arg, x, istart, nblock);
y = ei.shape(y);


% Gtomo2_table_block_back()
function x = Gtomo2_table_block_back(arg, y, istart, nblock)

ia = istart:nblock:arg.sg.na;
[y, eo] = embed_out(y, [arg.sg.nb length(ia)]);
x = Gtomo2_table_block_back_array(arg, y, istart, nblock);
x = eo.shape(x, arg.ig.mask, arg.ig.np);


%
% fatrix2 versions _array
%

% Gtomo2_table_block_forw_array()
function y = Gtomo2_table_block_forw_array(arg, x, istart, nblock)

ia = istart:nblock:arg.sg.na;
if nblock == 1 || size(arg.table,3) == 1
	table = arg.table;
else
	table = arg.table(:,:,ia);
end

y = jf_mex(arg.proj_str, arg.ptab_arg{:}, table, ...
	arg.mask8, arg.angles(ia), single(x));


% Gtomo2_table_block_back_array()
function x = Gtomo2_table_block_back_array(arg, y, istart, nblock)

ia = istart:nblock:arg.sg.na;
if nblock == 1 || size(arg.table,3) == 1
	table = arg.table;
else
	table = arg.table(:,:,ia);
end

x = jf_mex(arg.back_str, arg.ptab_arg{:}, table, ...
	arg.mask8, arg.angles(ia), single(y));
