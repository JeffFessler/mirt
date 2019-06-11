 function ftab = de_ftab_build(s_arg, varargin)
%function ftab = de_ftab_build(s_arg, [options])
%|
%| For multiple-kVp X-ray imaging, we must evaluate functions f_m(s_1, ..., s_L)
%| for m=1,..,M, where M is the number of kVp settings,
%| and L is the number of material components.
%| f_m(s1, s2) = -log( \int exp(- (m1(E) s1 + m2(E) s2)) dI_m(E) / I_m(E) )
%| (For "dual-energy" imaging, M=L=2.)
%|
%| This routine builds tables/models of f_m, its inverse, and its derivatives.
%| These are needed for (dual-energy) X-ray CT reconstruction.
%|
%| in:	(these all have sensible defaults, so try calling with no args)
%|	s_arg	cell	arguments for xray_read_spectra
%| option
%|	'mtype	cell	material types, e.g., {'soft', 'bone'}
%|	'ftype	char	fit type for de_ftab_fit()
%|	'sl	cell{LL} sample thicknesses for each of LL materials	
%|	's_n	[L,1]	number of samples of l'th material integrals
%|	's_max	[L,1]	maximum material density line integrals, units: g/cm^2
%|	'wt_fm	{M}	weighting for fitting fm
%|
%| out:
%|	ftab	struct
%|	methods:
%|		ftab.fit.fmfun(s1, s2, ...)
%|
%| Copyright 2001-04-27, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(s_arg, 'test'), de_ftab_build_test, return, end

ftab.show = false;
ftab.sl = {};
ftab.s_n = [];
ftab.s_min = [];
ftab.s_max = [];
ftab.wt_fm = {};
ftab.stop_after_fit = false;
ftab.mtype = {'soft', 'bone'}; % default materials
ftab.ftype = ''; % defer to de_ftab_fit() default
ftab = vararg_pair(ftab, varargin);

if isempty(s_arg)
%	s_arg = {'mono,60,100'}; % monoenergetic for testing
	s_arg = {'ps1'};
% fdir = sprintf('../fig,%s/', s_arg{1});
end

% material integral sampling:
% s_l is samples of the integral of the density of the lth material type.
[ftab.sl ftab.s_n ftab.s_min ftab.s_max] = ...
	de_ftab_build_sl(ftab.sl, ftab.s_n, ftab.s_min, ftab.s_max);

% Load X-ray spectra, etc.
ftab.xray = xray_read_spectra(s_arg{:}, 'show', ftab.show);
MM = size(ftab.xray.sp, 2);	% # of spectra
if ftab.show, prompt, end

% mac: mass attenuation coefficients
[ftab.mac LL] = de_ftab_build_mac(ftab.mtype, ...
	ftab.xray.en, ftab.xray.Ide, ftab.xray.eff, ftab.show);

% build tables of f_m(s1, s2)
% f_m looks very linear but it is not quite linear!
ftab.fm_fun = @de_ftab_build_fm;
ftab.fm = ftab.fm_fun(ftab, ftab.sl{:});
if ftab.show % figure showing f1, f2
	de_ftab_build_show_fm(ftab.fm, ftab.sl, 4)
prompt
end

% parametric fit to each fm to form a continuous function!
ftab.fit = de_ftab_fit(ftab.sl, ftab.fm, 'type', ftab.ftype, ...
	'mtype', ftab.mtype, 'wt', ftab.wt_fm);
%ftab.fit = de_ftab_fit(ftab.sl, ftab.fm, 'type', 'poly', 'order', 3, 'dc', 1);
if ftab.show % compare fit to sampled fm
	ftab.fit.show_fm(ftab.sl, ftab.fm)
	ftab.fit.show_err(ftab.sl, ftab.fm)
end

% for safefty, remove table, forcing user to use ftab.fit.fmfun!
% ftab = rmfield(ftab, 'fm');

if ftab.stop_after_fit, return, end

warning 'todo: not done after here'
return

%
% Build inverse of soft-tissue component
% for conventional "water only" beam-hardening correction.
%
if ~isvar('ftab.inv_water'), printm 'water inv'
	fw = ftab.feval(ftab, ftab.sl{1}, 0);	% mm=1 spectrum, water only
	fw = squeeze(fw(:,1,:));	% [s_n(1) M]
	ftab.inv_water_f = linspace(0, max(fw(:)), 51)';
	for mm=1:MM
		ftab.inv_water(:,mm) = interp1(fw(:,mm), ftab.sl{1}, ...
			ftab.inv_water_f, 'cubic', 'extrap');
	end
	if any(isnan(ftab.inv_water)), error 'nan', end

	ftab.inv_water_eval = 'interp1(ftab.inv_water_f, ftab.inv_water(:,m), f, ''cubic'', ''extrap'')';
	ftab.inv_water_eval = @(ftab,m,f) ftab.inv_water_eval(ftab, m, f);

	if ftab.show
		clf, subplot(221)
		plot(ftab.sl{1}, fw(:,1), '-', ...
			ftab.inv_water(:,1), ftab.inv_water_f, '.')
		axis tight, xlabel 's1', ylabel 'f1', title 'water fit m=1'
		subplot(222)
		plot(ftab.sl{1}, fw(:,2), '-', ...
			ftab.inv_water(:,2), ftab.inv_water_f, '.')
		axis tight, xlabel 's1', ylabel 'f2', title 'water fit m=2'

		wat(:,1) = ftab.inv_water_eval(ftab, 1, fw(:,1));
		wat(:,2) = ftab.inv_water_eval(ftab, 2, fw(:,2));
		err = wat - [ftab.sl{1} ftab.sl{1}];

		subplot(223)
		plot(ftab.sl{1}, err(:,1), '.-')
		axis tight, xlabel s1, ylabel err1, title 'Water inverse error'
		subplot(224)
		plot(ftab.sl{1}, err(:,2), '.-')
		axis tight, xlabel s1, ylabel err2, title 'Water inverse error'

		prompt
	end
	clear fw t2 wat err
end


%
% Jacobian transformation (not essential; for exploration only)
%
if ~isvar('ftab.T'), printm 'jacobian f->fs approx s!'
	ftab.T = pinv(mass.bar);	% [LL,MM] transformation, fs = T f

	ftmp = ftab.feval(ftab, sll{:});
	if ftab.show
		zz = de_ftab_xform(ftab, ftmp);
		z1 = reshape(zz(:,:,1), ftab.s_n);
		z2 = reshape(zz(:,:,2), ftab.s_n);
	end

	if 0 % array of dots showing effect of transformation
		clf, plot(z1, z2, '.')
		axis([0 max(ftab.sl{1}) 0 max(ftab.sl{2})])
		xlabel 's_1 [cm^2/g]', ylabel 's_2 [cm^2/g]'
		xlabel 'f^*_1(s_1,s_2)', ylabel 'f^*_2(s_1,s_2)'
		title 'Effect of Linearizing Transformation'
	return
	end

%	ir_savefig('c', fdir, 'fig_f_linearized')
	clear ftmp zz % can't clear z1,z2 yet!

	%
	% nice figure showing inverse and departure from linearity thereof
	%
	if ftab.show
		t.i1 = 1:4:length(ftab.sl{1});
		t.i2 = 1:4:length(ftab.sl{2});
		t.z1 = z1(t.i1, t.i2);
		t.z2 = z2(t.i1, t.i2);
		t.s1 = sll{1}(t.i1, t.i2);
		t.s2 = sll{2}(t.i1, t.i2);
		clf
%		plot3(t.z1, t.z2, t.s1-0*z1, '.')
		subplot(221), mesh(t.z1, t.z2, t.s1-0*t.z1), title 's_1'
		axis([0 ftab.s_max(1) 0 ftab.s_max(2) 0 ftab.s_max(1)])
		xtick([]), ytick([]), zwhite, % xlabel f^*_1, ylabel f^*_2

		subplot(222), mesh(t.z1, t.z2, t.s2-0*t.z2), title 's_2'
		axis([0 ftab.s_max(1) 0 ftab.s_max(2) 0 ftab.s_max(2)])
		xtick([]), ytick([]), ztick, zwhite
		% xlabel f^*_1, ylabel f^*_2

		subplot(223), mesh(t.z1, t.z2, t.s1-1*t.z1), title 's_1 - f^*_1'
		axis([0 ftab.s_max(1) 0 ftab.s_max(2) -3 12])
		xtick, ytick, ztick, zwhite
		% xlabel f^*_1, ylabel f^*_2

		subplot(224), mesh(t.z1, t.z2, t.s2-1*t.z2), title 's_2 - f^*_2'
		axis([0 ftab.s_max(1) 0 ftab.s_max(2) -3 12])
		xtick, ytick, ztick, zwhite, xlabel f^*_1, ylabel f^*_2

%		ir_savefig('c', fdir, 'fig_fs_inverse')
		clear t

		prompt
	end
end


%
% Build polynomial approximation to shat = (T F)^{-1}(z).
% Use weighted fitting so that soft tissue is fit well!
%
if ~isvar('ftab.inv.eval'), printm 'inv: f -> s'
	if 0	% this polynomial way seems insufficiently accurate!?
		ftab.inv.basis_func = ir_poly2_fun(3);
	%	wt = sll{2}.^3;	% weighting
		wt = 1 ./ (0*sll{2} + 1);	% weighting
		tmp = ftab.inv.basis_func(z1(:), z2(:));
		ftab.inv.nbasis = length(ftab.inv.basis_func(0,0));
		ftab.inv.coef = zeros(ftab.inv.nbasis,2);
		for ll=1:LL
			ftab.inv.coef(:,ll) = (diag_sp(wt(:)) * tmp) \ ...
				(wt(:) .* sll{ll}(:));
		end
	end

	ftab.inv.eval = @(ftab, fhat) de_ftab_invert(ftab,fhat);

	% test inverse with a different set of samples to see errors
	if 1
		s1 = linspace(0, ftab.s_max(1)+1, length(ftab.sl{1})+8);
		s2 = linspace(0, ftab.s_max(2)+1, length(ftab.sl{2})+8);
		ss = zeros(length(s1), length(s2), 2);
		[ss1 ss2] = ndgrid(s1, s2);
		ftmp = ftab.feval(ftab, ss1, ss2);
		stmp = ftab.inv.eval(ftab,ftmp);

		err = stmp - stackup(ss1,ss2);
		err = abs(err);
		printm('worst inverse error: %g of %g', max(err(:)), max(s1(:)))
	end

	if ftab.show
		clf, subplot(221), mesh(s1, s2, ss1'), title 's_1'
		axis tight, xtick, ytick, zwhite, grid
		hold on, plot3(ss1, ss2, stmp(:,:,1), 'y.'), hold off
		xlabel s_1, ylabel s_2

		subplot(222), mesh(s1, s2, ss2'), title 's_2'
		axis tight, xtick, ytick, zwhite, grid
		hold on, plot3(ss1, ss2, stmp(:,:,2), 'y.'), hold off

		subplot(223), mesh(s1, s2, err(:,:,1)')
		colormap hsv, caxis([0 max(err(:))]), cbar
		axis tight, xtick, ytick, zwhite, grid
		subplot(224), mesh(s1, s2, err(:,:,2)')
		axis tight, xtick, ytick, zwhite, grid
		xlabel 's_1', ylabel 's_2', title 'error 2'

		prompt
	end
	clear t z1 z2 tmp wt
end


if 0
	name = ['table/' s_arg{1}];
	ans = sprintf('save table to "%s"? ', name);
	ans = input(ans, 's');
	if isempty(ans) | ~streq(ans, 'y'), return, end
	save(name, 'ftab')
	printm('table saved to "%s"', name)
end

end % de_ftab_build()


%
% de_ftab_build_sl()
%
function [sl, s_n, s_min, s_max] = de_ftab_build_sl(sl, s_n, s_min, s_max)

% number of samples of the material integrals "s"
if isempty(s_n)
	if isempty(sl)
		s_n = [45 43];
	else
		for ll=1:length(sl)
			s_n(1,ll) = length(sl{ll});
		end
	end
end

% minimum material "integrals"
if isempty(s_min)
	if isempty(sl)
		s_min = [0 0];
	else
		for ll=1:length(sl)
			s_min(1,ll) = min(sl{ll});
		end
	end
end

% maximum material "integrals"
if isempty(s_max)
	if isempty(sl)
		% soft max: 50cm * 1g/cc
		% bone max: 15cm * 2g/cc (for now)
		s_max = [50 30];
	else
		for ll=1:length(sl)
			s_max(1,ll) = max(sl{ll});
		end
	end
end

if isempty(sl)
	for ll=1:length(s_n)
		sl{ll} = linspace(s_min(ll), s_max(ll), s_n(ll))';
	end
end

end % de_ftab_build_sl()


%
% de_ftab_build_mac()
%
function [mac LL] = de_ftab_build_mac(type, en, Ide, xeff, show)
if streq(type{1}, 'pca_given')
	mac.type = {'pca1', 'pca2'}; % hardwired to first two components
	mac.mac = type{2}.basis; % [ne,2]

elseif streq(type{1}, 'pca_do')
	mac.type = {'pca1', 'pca2'}; % hardwired to first two components
	tmp = type{2}; % arguments for de_component2()
	tmp = de_component2('kev', en, tmp{:});
	mac.mac = tmp.basis; % [ne,2]

else
	mac.type = type;
	% interpolate mass atten coef to source energy sampling
	mac.mac = xray_read_atten(type, en); % [ne,L]
end

mac.bar = diag(1 ./ sum(Ide)) * (Ide' * mac.mac);

LL = length(type);

% examine condition number
if 1
	printm 'bmassml = '
	disp(mac.bar)
	% inv(mass.bar)
%	printm('condition = %g', cond(mass.bar))
	pr cond(mac.bar)
%	pr inv(mac.bar' * mac.bar)
end

if show
	clf
	plot(	en, mac.mac(:,1), 'c:', ...
		en, mac.mac(:,2), 'y:', ...
		xeff(1), mac.bar(1,:), 'g+', ...
		xeff(2), mac.bar(2,:), 'r^')
	axisy(0.1, max(0.7, 1.1*max(mac.bar(:))))
	legend(mac.type{:})
prompt
end

end % de_ftab_build_mac()


%
% de_ftab_build_fm()
% build tables of f_m(s1, s2, ...)
%
function fm = de_ftab_build_fm(ftab, varargin)
if length(varargin) == 1 && iscell(varargin{1})
	sll = varargin{1};
else
	sll = ndgrid_jf('cell', varargin{:});
end
fm = de_ftab_fm(sll, ftab.mac.mac, ftab.xray.Ide);

end % de_ftab_build_fm()


%
% de_ftab_build_show_fm()
%
function de_ftab_build_show_fm(fm, sl, down)
fmax = 1.01 * max(fm(:));
fmax = ceil(fmax);
s_max(1) = max(sl{1}(:));
s_max(2) = max(sl{2}(:));

i1 = 1:down:length(sl{1});
i2 = 1:down:length(sl{2});
s1 = sl{1}(i1);
s2 = sl{2}(i2);
f1 = fm(i1,i2,1);
f2 = fm(i1,i2,2);

clf
subplot(221), mesh(s1, s2, f1')
colormap hsv, caxis([0 fmax]), cbar
axis([0 s_max(1) 0 s_max(2) 0 fmax])
xtick, ytick, ztick, zwhite, xlabel s_1, ylabel s_2, title 'f_1(s)'
text(-60, 10, '[cm^2/g]')

subplot(222), mesh(s1, s2, f2')
colormap hsv, caxis([0 fmax]), cbar
axis([0 s_max(1) 0 s_max(2) 0 fmax])
xtick, ytick, ztick, zwhite, xlabel s_1, ylabel s_2, title 'f_2(s)'

% ir_savefig('c', fdir, 'fig_f1_f2')
end % de_ftab_build_show_fm()


%
% de_ftab_build_test
%
function de_ftab_build_test
ftab = de_ftab_build({'ps1'}) % , 'mtype', {'pca_do', {}})
end % de_ftab_build_test
