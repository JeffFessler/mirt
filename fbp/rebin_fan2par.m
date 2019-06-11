  function psino = rebin_fan2par(fsino, sf, sp, varargin)
%|function psino = rebin_fan2par(fsino, sf, sp, varargin)
%|
%| Rebin fan-beam sinogram into parallel-beam (or mojette) sinogram.
%| (Also useful for parallel-to-mojette rebinning.)
%|
%| in
%|	fsino	[ns nbeta (L)]	fan-beam sinogram (or possibly parallel)
%|	sf			sino_geom() for fan-beam input
%|	sp			sino_geom() for parallel-beam output
%|
%| options
%|	'ob'			set to 1 to create (fatrix2 or Fatrix) object
%|	'class'			'fatrix2' (default) or 'Fatrix' (obsolete)
%|	's_interp'	default: {'order', 3, 'ending', 'zero'}
%|	'beta_interp'	default: {'order', 3, 'ending', 'periodic'}
%|
%| out
%|	psino	[nr nphi (L)]	parallel-beam sinogram (or possibly mojette)
%|
%| Copyright 2005-12-10, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fsino, 'test'), rebin_fan2par_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.ob = false;
arg.class = 'fatrix2';
arg.s_interp = {'order', 3, 'ending', 'zero'};
arg.beta_interp = {'order', 3, 'ending', 'periodic'};
arg.na_min_fail = 180; % warn if fewer than this many views
arg = vararg_pair(arg, varargin);

if isempty(sf.nb), sf.nb = size(fsino,1); end
if isempty(sf.na), sf.na = size(fsino,2); end
arg.dimi = sf.dim;
arg.dimo = sp.dim;

if (sf.na < arg.na_min_fail) % angular interpolation needs decent sampling!
	fail('na=%d too few; need %d', sf.na, arg.na_min_fail)
end

if streq(sp.type, 'par')
	is_mojette = 0;
	dr = sp.dr;
elseif streq(sp.type, 'moj')
	is_mojette = 1;
	dr = sp.dx; % trick:
else
	error 'output sino must be "par" or "moj"'
end

arg.sf = sf;

if streq(sf.type, 'fan') % fan->(par|moj)
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_fan2par_setup(...
		sf.nb, sf.ds, sf.offset_s, ...
		sf.na, sf.orbit_start, sf.orbit, ...
		sf.dso, sf.dsd, sf.dfs, ...
		sp.nb, dr, sp.offset_r, ...
		sp.na, sp.orbit_start, sp.orbit, ...
		is_mojette, arg.s_interp, arg.beta_interp, arg.class);

elseif streq(sf.type, 'par') % par->(moj|par)
	[arg.r_ob arg.phi_ob arg.flag180] = rebin_fan2par_setup(...
		sf.nb, sf.dr, sf.offset_r, ...
		sf.na, sf.orbit_start, sf.orbit, ...
		inf, inf, 0, ...
		sp.nb, dr, sp.offset_r, ...
		sp.na, sp.orbit_start, sp.orbit, ...
		is_mojette, arg.s_interp, arg.beta_interp, arg.class);

else
	fail([sp.type ' to ' sf.type ' not done'])
end

if arg.ob
	psino = rebin_fan2par_ob(arg);
else
	arg.nb = sp.nb;
	arg.na = sp.na;
	psino = rebin_fan2par_arg_many(arg, fsino);
end


% rebin_fan2par_ob()
function ob = rebin_fan2par_ob(arg)

switch arg.class
case 'Fatrix'
	dim = [prod(arg.dimo) prod(arg.dimi)];
	ob = Fatrix(dim, arg, 'caller', 'rebin_fan2par', ...
        	'forw', @rebin_fan2par_arg_Fatrix, 'back', @error);
case 'fatrix2'
	ob = fatrix2('arg', arg, 'does_many', 1, ...
		'odim', arg.dimo, 'idim', arg.dimi, ...
        	'forw', @rebin_fan2par_arg, 'back', @error);
otherwise
	fail 'bad class'
end


% rebin_fan2par_arg()
function sino = rebin_fan2par_arg(arg, sino)

sino = (arg.phi_ob * sino.').';

if arg.flag180
	sino = rebin_fan2par_inlace(arg, sino);
end

sino = arg.r_ob * sino;


% rebin_fan2par_arg_many()
function psino = rebin_fan2par_arg_many(arg, fsino)
tmp = size(fsino);
nt = prod(tmp(3:end));
fsino = reshape(fsino, [tmp(1) tmp(2) nt]);
psino = zeros(arg.nb, arg.na, nt, class(fsino));
for it=1:nt
	psino(:,:,it) = rebin_fan2par_arg(arg, fsino(:,:,it));
end
psino = reshape(psino, [arg.nb arg.na tmp(3:end)]);


% rebin_fan2par_arg_Fatrix()
function sino = rebin_fan2par_arg_Fatrix(arg, sino)

if size(sino,1) == prod(arg.dimi)
	sino = reshape(sino, arg.dimi(1), arg.dimi(2), []);
	flag_column = 1;
else
	flag_column = 0;
end

sino = rebin_fan2par_arg(arg, sino);

if flag_column
	sino = reshape(sino, prod(arg.dimo), []);
end


% rebin_fan2par_inlace()
% interlace opposing views (for quarter-detector offset)
function sino = rebin_fan2par_inlace(arg, sino)

if ndims(sino) > 2, error 'multisino not done', end
ns = arg.dimi(1);
nphi = arg.dimo(2);
t1 = sino(:,1:nphi); % [ns nphi]
t2 = flipdim(sino(:,nphi + [1:nphi]), 1); % [ns nphi]
if arg.sf.offset == 1.25 % trick
	t1 = [t1 ; t2([ end-1 end],:)]; % [ns+2 nphi]
	t2 = [t1([1 2],:) ; t2]; % [ns+2 nphi]
	ns = ns + 2;
elseif arg.sf.offset ~= 0.25
	error 'bug'
end
sino = reshape([t1(:)'; t2(:)'], 2*ns, nphi, []);


% rebin_fan2par_setup()
function [r_ob phi_ob flag180] = rebin_fan2par_setup(...
	ns, ds, offset_s, ...
	nbeta, beta_start, beta_orbit, ...
	dso, dsd, dfs, ...
	nr, dr, offset_r, ...
	nphi, phi_start, phi_orbit, ...
	is_mojette, ... % 1 if output sinogram is to be mojette
	s_interp, beta_interp, ob_class)

if dfs, error 'flat fan not done', end

if phi_orbit == 180 && beta_orbit == 360
	flag180 = 1; % trick: handle 180 -> 360
else
	flag180 = 0;
end

phi_start = deg2rad(phi_start);
phi_orbit = deg2rad(phi_orbit);
beta_start = deg2rad(beta_start);
beta_orbit = deg2rad(beta_orbit);
phi = phi_start + phi_orbit * [0:nphi-1]' / nphi;
%bet = beta_start + beta_orbit * [0:nbeta-1]' / nbeta;

% angular interpolator
ws = (ns-1)/2 + offset_s;

if isinf(dsd) % parallel - no need for angular interp if orbits match
	if phi_start == beta_start && phi_orbit == beta_orbit && nphi == nbeta
		phi_ob = 1;
	else
		error 'parallel beam angular interpolation not done'
	end

else % fan

	if flag180 % expand desired phi's from [x,x+180] to [x,x+360]
		phi2 = phi_start + 2 * phi_orbit * [0:(2*nphi-1)]' / (2*nphi);
	else
		if phi_orbit ~= deg2rad(360) || beta_orbit ~= deg2rad(360)
			error 'todo: only 360 done - ask jeff'
		end
		phi2 = phi;
	end

	s = ([0:ns-1]' - ws) * ds;

	bet = outer_sum(phi2, -s / dsd); % beta = phi - gam
	bet_int = nbeta / beta_orbit * (bet - beta_start);

	phi_ob = bspline_1d_interp(zeros(nbeta,ns), bet_int, ...
		beta_interp{:}, 'ob', 1, 'class', ob_class);
end

% radial interpolator

wr = (nr-1)/2 + offset_r;
if is_mojette
	dr = dr * max(abs(cos(phi)), abs(sin(phi)))';
end
r = ([0:nr-1]' - wr) * dr; % trick: [nr 1] or [nr nphi]
if isinf(dsd)
	s = r;
else
	s = dsd * asin(r / dso);
end

if flag180 % trick: due to interlacing, the effective "s" sampling changes
	if offset_s ~= 0.25 && offset_s ~= 1.25
		error 'only 0.25 and 1.25 implemented now; but generalizable'
	end
	ns = ns * 2 + 4 * (offset_s - 0.25);
	offset_s = 0;
	ws = (ns-1)/2 + offset_s;
	ds = ds / 2;
end
s_int = s / ds + ws;
r_ob = bspline_1d_interp(nan(ns, nphi), s_int, s_interp{:}, ...
	'ob', 1, 'class', ob_class);


% rebin_fan2par_test
% test both rebin_fan2par() and par2fan_rebin()
function rebin_fan2par_test

if 1
	rebin_fan2par_plot
	prompt
end

down = 4;
dx = down/2;

% todo: fix dr to use formula in 3.9.1 !
gp = sino_geom('par', 'nb', 1096/down, 'na', 800/down, ...
	'strip_width', 'd', ...
	'dr', 0.5*down, 'offset_r', 0, 'orbit', 180);
gm = sino_geom('moj', 'nb', gp.nb, 'na', gp.na, ...
	'strip_width', 'd', ...
	'dx', dx, 'offset_r', 0, 'orbit', 180);
gf = sino_geom('fan', 'ns', 888/down, 'nbeta', 984/down, ...
	'strip_width', 'd', ...
	'dsd', 949, 'dod', 408, ...
	'ds', down, 'offset_s', 1.25); % quarter detector

% analytical sinograms
ell = [0 20 180 150 0 1];
oversample = 4;
par = ellipse_sino(gp, ell, 'oversample', oversample);
moj = ellipse_sino(gm, ell, 'oversample', oversample);
fan = ellipse_sino(gf, ell, 'oversample', oversample);

% sanity check for par->par
if 0
	g0 = sino_geom('fan', 'nb', gp.nb, 'na', gp.na, ...
		'dsd', inf, 'dod', 0, 'orbit', gp.orbit, ...
		'ds', gp.dr, 'offset_s', gp.offset_r);
	par2par = rebin_par2fan(par, gp, g0);
	max_percent_diff(par, par2par)
	par2par = rebin_fan2par(par, g0, gp);
	max_percent_diff(par, par2par)
return
end

if 1 % test rebin (normal)
	cpu etic
	fan2par = rebin_fan2par(fan, gf, gp);
	cpu etoc 'fan2par time'

	cpu etic
	par2fan = rebin_par2fan(par, gp, gf);
	cpu etoc 'par2fan time'

	max_percent_diff(par, fan2par)
	max_percent_diff(fan, par2fan)

	im plc 2 3
	im(1, par, 'par'), cbar
	im(4, fan, 'fan'), cbar
	im(2, fan2par, 'fan2par'), cbar
	im(5, par2fan, 'par2fan'), cbar
	im(3, fan2par-par, 'error'), cbar
	im(6, par2fan-fan, 'error'), cbar

	if 1 % check object
		classes = {'Fatrix', 'fatrix2'};
		for ic = 1:numel(classes)
			tmp = rebin_fan2par(nan(size(fan)), gf, gp, ...
				'ob', 1, 'class', classes{ic});
			tmp = tmp * fan;
			jf_equal(tmp, fan2par)

			tmp = rebin_par2fan(nan(size(par)), gp, gf, ...
				'ob', 1, 'class', classes{ic});
			tmp = tmp * par;
			jf_equal(tmp, par2fan)
		end
	end
prompt
end

% test rebin (mojette)
cpu etic
fan2moj = rebin_fan2par(fan, gf, gm);
cpu etoc 'fan2moj time'

cpu etic
moj2fan = rebin_par2fan(moj, gm, gf);
cpu etoc 'moj2fan time'

max_percent_diff(fan, moj2fan)
max_percent_diff(moj, fan2moj)

im(1, moj, 'moj'), cbar
im(4, fan, 'fan'), cbar
im(2, fan2moj, 'fan2moj'), cbar
im(5, moj2fan, 'moj2fan'), cbar
im(3, fan2moj-moj, 'error'), cbar
im(6, moj2fan-fan, 'error'), cbar


% rebin_fan2par_plot
function rebin_fan2par_plot

down = 8;
sg = sino_geom('ge1', 'down', down);
dsd = sg.dsd;

xds = sg.xds;
yds = sg.yds;
d0 = [xds yds];
c0 = [0 -sg.dod];
clf
s0 = [0 sg.dso];
plot(0, 0, '.', s0(1), s0(2), 'o', d0(:,1), d0(:,2), '.-')
axis square
axis equal
hold on
plot([s0(1) c0(1)], [s0(2) c0(2)], ':')
hold off

blist = linspace(-1,1,7) * pi/8;
nb = numel(blist);
for ii=1:nb
	bet = pi/8;
	bet = blist(ii);
	rot = [cos(bet) sin(bet); -sin(bet) cos(bet)];

	cr = c0 * rot;
	dr = d0 * rot;
	sr = s0 * rot;

	hold on
	h = plot(sr(1), sr(2), 'o', dr(:,1), dr(:,2), '.-');
%	get(h(1))
	col = [0 0 ii/nb];
	col = colormap(hsv(nb+3));
	col = col(ii,:);
	set(h(1), 'color', col)
	set(h(2), 'color', col)
%keyboard
	h = plot([sr(1) cr(1)], [sr(2) cr(2)], ':');
	set(h, 'color', col)
	h = plot([1 1]*sr(1), sr(2)-[0 dsd], '--');
	set(h, 'color', col)
	hold off
end
xlabelf 'x'
ylabelf 'y'
% ir_savefig cw fig_rebin_fan2par
