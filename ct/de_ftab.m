 function ftab = de_ftab(xrs, mas, varargin)
%function ftab = de_ftab(xrs, mas, [options])
%|
%| For multiple-kVp X-ray imaging, we must evaluate functions f_m(s_1, ..., s_L)
%| for m=1,..,M, where M is the number of kVp settings,
%| and L is the number of material components.
%| f_m(s1, s2) = -log( \int exp(- (m1(E) s1 + m2(E) s2)) dI_m(E) / I_m(E) )
%| (For "dual-energy" imaging, M = L = 2.)
%|
%| This routine builds tables/models of f_m, its inverse, and its derivatives.
%| These are needed for (dual-energy) X-ray CT reconstruction.
%|
%| in
%|	xrs	strum	X-ray spectra; see xray_read_spectra.m
%|	mas	strum	mass attenuation coefficients; see xray_read_mac.m
%|
%| option
%|	'sls'	strum	see de_ftab_sls.m
%|	'ftype'	char	fit type for de_ftab_fit() (default: '')
%|	'fit_args' {}	options for de_ftab_fit() (default: {})
%|			e.g.: {'wt', num2cell(ones(M,1))} weighting for fitting
%|	'ctype'	char	curv type for de_ftab_curv() (default: '')
%|	'show'	0|1	visualize
%|
%| out
%|	ftab	strum
%|		data:
%|			.xrs,.mas,.sls,.mac	based on input
%|			.fit	strum with fm approximation methods
%|				see de_ftab_fit.m
%|			.inv1	strum to invert 1 component BH, e.g., water
%|				see de_ftab_inv1.m
%|			.inv2	strum to invert BH by polynomials,
%|				see de_ftab_inv1.m
%|		methods:
%|			.fm_fun(s1, s2, ...) (probably not needed by user)
%|			.plot_fm()
%|			.plot_mac()
%|			.plot_jac()
%|			.plot_inv()
%|
%| Copyright 2008-6-15, Jeff Fessler, University of Michigan

if nargin == 1 && streq(xrs, 'test'), de_ftab_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

ftab.show = false;
ftab.sls = struct([]);
ftab.ftype = ''; % defer to de_ftab_fit() default
ftab.fit_args = {}; % defer to de_ftab_fit() default(s)
ftab.ctype = ''; % defer to de_ftab_curv() default ('newt' is fastest)
ftab = vararg_pair(ftab, varargin);

% X-ray spectra
if isempty(xrs)
%	xrs = {'mono,60,100'}; % monoenergetic for testing
	xrs = {'ps1'};
end
if iscell(xrs)
	xrs = xray_read_spectra(xrs{:});
end
if ftab.show
	clf, xrs.plot
prompt
end
ftab.xrs = xrs;

ftab.MM = size(xrs.sp, 2); % # of spectra

% mas: mass attenuation coefficient (strum), at tabulated energies
if isempty(mas)
	mas = {'soft', 'bone'};	% default materials
end
if iscell(mas)
	mas = xray_read_mac(mas);
end
if ftab.show
	clf, mas.plot('kev', xrs.en)
prompt
end
ftab.mas = mas;

% mac: mass attenuation coefficient, at spectrum energies
ftab.mac = xray_make_mac(xrs, mas);
ftab.LL = ncol(ftab.mac.mac);
if ftab.show
	pr ftab.mac.bar % 'bmassml'
%	pr cond(ftab.mac.bar) % examine condition number
	pr sqrt(cond(ftab.mac.bar' * ftab.mac.bar)) % examine condition number
	de_ftab_plot_mac(ftab);
prompt
end

% material integral sampling:
% s_l is samples of the integral of the density of the lth material type.
if isempty(ftab.sls)
	ftab.sls = de_ftab_sls;
end

% build tables of f_m(s1, s2)
fm_fun = @de_ftab_make_fm;
ftab.fm = fm_fun(ftab);
if ftab.show % figure showing f1, f2.  it looks linear, but it not quite!
	de_ftab_plot_fm(ftab, 'down', 4);
prompt
end

% parametric fit to each fm to form a continuous function!
printm 'fit'
if streq(ftab.ftype, 'exp3')
	ftab.fit = de_ftab_fit_exp3(ftab.xrs, ftab.mac.mac, ftab.fit_args{:});
else
	ftab.fit = de_ftab_fit(ftab.sls.sl, ftab.fm, 'type', ftab.ftype, ...
		'mtype', mas.type, ...
		'kev', xrs.en, 'mac', ftab.mac.mac, ... % needed for exp
		'macbar', ftab.mac.bar, ftab.fit_args{:});
% todo: may be bad (and slow) to use finely sampled energies (kev) for exp fit
%		'type', 'exp');
%ftab.fit = de_ftab_fit(ftab.sl, ftab.fm, 'type', 'poly', 'order', 3, 'dc', 1);
end

ftab.fit = de_ftab_curv(ftab.fit, 'ctype', ftab.ctype);

% Jacobian transformation
ftab.T = pinv(ftab.mac.bar); % [L M] transformation, fs = T f approx s !

% Build 1D inverse of 1st material component (usually water)
% for conventional "water only" beam-hardening correction.
ftab.inv1 = de_ftab_inv1(ftab.fit, ftab.sls.sl{1});

% build polynomial inverse approximation
ftab.inv2 = de_ftab_inv2(ftab.fit, ftab.sls.sl, 'T', ftab.T); % todo: options

%ftab.inv.eval = @(ftab,fhat) de_ftab_invert(ftab,fhat);
%ftab.inv = de_ftab_inv_setup(ftab);

meth = {'fm_fun', fm_fun, '(cell: sl)';
	'show_fm', @de_ftab_plot_fm, '()'; % backward compat
	'plot_fm', @de_ftab_plot_fm, '()';
	'plot_mac', @de_ftab_plot_mac, '()';
	'plot_jac', @de_ftab_plot_jac, '()';
	};

ftab = strum(ftab, meth);

ftab = de_ftab_inv_setup(ftab);

if ftab.show % compare fit to sampled fm, etc.
	ftab.inv1.plot(ftab.fit);
	if ~isempty(ftab.inv2)
		ftab.inv2.plot(ftab.fit, ftab.sls.sl);
	end
%	ftab.fit.show_fm(ftab.sls.sl, ftab.fm)
%	ftab.fit.show_err(ftab.sls.sl, ftab.fm)
	ftab.plot_jac;
	ftab.plot_inv;
%keyboard
end

% for safefty, remove table, forcing user to use ftab.fm_fun or ftab.fit.fmfun !
% ftab = rmfield(ftab, 'fm');

end % de_ftab()


%
% de_ftab_inv_setup()
% add methods related to inverse to ftab
%
function ftab = de_ftab_inv_setup(ftab)

printm 'inv: f -> s'

arg.inv = struct; % place holder
fun = @(ftab, fhat) de_ftab_invert(ftab, fhat);
meth = {'inv_fun', fun, '(fhat)';
	'plot_inv', @de_ftab_plot_inv, '()'};
ftab = strum(arg, meth, 'base', ftab);

end % de_ftab_inv_setup()


%
% de_ftab_plot_inv()
% plot inverse, testing with a different set of samples to see errors
%
function dummy = de_ftab_plot_inv(ftab)
dummy = [];
LL = ftab.LL;
MM = ftab.MM;
if LL ~= 2, warn 'L=2 done only', return, end
if MM < LL, warn('plot_inv with MM=%d < LL=%d skipped', MM, LL), return, end

if 1
	sls = ftab.sls;
	s1 = linspace(0, sls.max(1)+1, length(sls.sl{1}) + 8);
	s2 = linspace(0, sls.max(2)+1, length(sls.sl{2}) + 8);
	ss = ndgrid_jf('mat', s1, s2);
	ftmp = ftab.fit.fmfun(ss);
	stmp = ftab.inv_fun(ftmp);

	err = stmp - ss;
	err = abs(err);
	for ll=1:LL
		printm('worst inverse error l=%d: %g of %g', ll, ...
			max(col(stackpick(err,ll))), max(col(stackpick(ss,ll))))
	end
end

if im && usejava('jvm')
	ss1 = ss(:,:,1);
	ss2 = ss(:,:,2);
	clf, pl = @(n) subplot(220 + n);
	pl(1), mesh(s1, s2, ss1')
	axis tight, xtick, ytick, zwhite, grid
	hold on, plot3(ss1, ss2, stmp(:,:,1), 'y.'), hold off
	xlabel 's1', ylabel 's2', title 's1'

	pl(2), mesh(s1, s2, ss2')
	axis tight, xtick, ytick, zwhite, grid
	hold on, plot3(ss1, ss2, stmp(:,:,2), 'y.'), hold off
	xlabel 's1', ylabel 's2', title 's2'

	pl(3), mesh(s1, s2, err(:,:,1)')
	colormap hsv, caxis([0 max(err(:))]), cbar
	axis tight, xtick, ytick, zwhite, grid
	xlabel 's1', ylabel 's2', title 'error 1'

	pl(4), mesh(s1, s2, err(:,:,2)')
	axis tight, xtick, ytick, zwhite, grid
	xlabel 's1', ylabel 's2', title 'error 2'
prompt
end

end % de_ftab_plot_inv()



%
% de_ftab_plot_mac()
%
function dummy = de_ftab_plot_mac(ftab, varargin)
dummy = [];
xrs = ftab.xrs;
mas = ftab.mas;
mac = ftab.mac;

ltype = {'c:', 'y:', 'g:', 'r:'};
mtype = {'g+', 'r^', 'm*', 'yo'};
arg = {};
for ll=1:ftab.LL
	arg = {arg{:}, xrs.en, mac.mac(:,ll), ltype{ll}};
end
for mm=1:ftab.MM
	arg = {arg{:}, xrs.eff(mm), mac.bar(mm,:), mtype{mm}};
end

plot(arg{:})
axisy(0.1, max(0.7, 1.1*max(mac.bar(:))))
xlabel 'E', ylabel 'mac(E)'
mtype = mas.type;
legend(mtype{:})

prompt
end % de_ftab_plot_mac()


%
% de_ftab_make_fm()
% build tables of f_m(s1, s2, ...)
%
function fm = de_ftab_make_fm(ftab, varargin)
if length(varargin) == 0
	varargin = ftab.sls.sl; % default
end
if length(varargin) == 1 && iscell(varargin{1})
	sll = varargin{1};
else
	sll = ndgrid_jf('cell', varargin{:});
end
fm = de_ftab_fm(sll, ftab.mac.mac, ftab.xrs.Ide);

end % de_ftab_make_fm()


%
% de_ftab_plot_fm()
%
function dummy = de_ftab_plot_fm(ftab, varargin)
dummy = [];
arg.down = 1;
arg = vararg_pair(arg, varargin);

fm = ftab.fm;
sl = ftab.sls.sl;

if ftab.LL == 1
	plot(sl{1}, fm, '-', sl{1}, sl{1} * ftab.mac.bar(:)', ':')
	xlabel 's1'
	ylabel 'f1(s1)'
	name = ftab.xrs.name;
	legend(name{:}, 2)

elseif ftab.LL == 2

	fmax = 1.01 * max(fm(:));
	fmax = ceil(fmax);
	s_max(1) = max(sl{1}(:));
	s_max(2) = max(sl{2}(:));

	i1 = 1:arg.down:length(sl{1});
	i2 = 1:arg.down:length(sl{2});
	s1 = sl{1}(i1);
	s2 = sl{2}(i2);

	pl = @(mm) subplot(1,ftab.MM,mm);

	clf
	for mm=1:ftab.MM
		f1 = fm(i1,i2,mm);
		pl(mm)
		plot(s1, f1', '.-')
		axisx([0 s_max(1)]), xlabel 's_1'
		titlef('f_%d(s)', mm)
		text(-60, 10, '[cm^2/g]')
	end

	if usejava('jvm')
		prompt
	else
		warn 'need jvm for plot_fm mesh'
	return
	end

	clf
	for mm=1:ftab.MM
		f1 = fm(i1,i2,mm);
		pl(mm)
		mesh(s1, s2, f1')
		colormap hsv, caxis([0 fmax]), cbar
		axis([0 s_max(1) 0 s_max(2) 0 fmax])
		xtick, ytick, ztick, zwhite, xlabel s_1, ylabel s_2
		titlef('f_%d(s)', mm)
		text(-60, 10, '[cm^2/g]')
	end
end

end % de_ftab_plot_fm()


%
% de_ftab_test
%
function de_ftab_test
xrs = xray_read_spectra('ps1'); % M=2
%xrs = xray_read_spectra('poly1,60,90,120'); % M=3 test
%xrs = xray_read_spectra('poly1,160', 'filters', {{'copper', 0.1}}); % M < L test
xrs.plot
prompt
mas = xray_read_mac({'water', 'bone'});
if 0
	mas = xray_read_mac({'soft', 'iodine'});
	tmp = mas.mac_raw;
	tmp{2}(:) = tmp{2}(:) / 10; % todo: explore scaling to improve cond #
	mas.mac_raw = tmp;
end
ftab = de_ftab(xrs, mas, 'show', im, ...
	 'ftype','exp', 'fit_args', {'kev', 10:5:160, 'mac', []});
if im
	ftab.plot_fm;
	prompt
	ftab.inv1.plot(ftab.fit);
	if ~isempty(ftab.inv2)
		ftab.inv2.plot(ftab.fit, ftab.sls.sl);
	end
	clf, ftab.plot_mac;
	ftab.plot_jac;
	ftab.plot_inv;
end
end % de_ftab_test
