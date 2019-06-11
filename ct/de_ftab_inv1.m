 function inv1 = de_ftab_inv1(fit, s1, varargin)
%function inv1 = de_ftab_inv1(fit, s1, [options])
%|
%| Build object that does 1D inverse of BH function for 1st material component,
%| (usually water), for conventional "water only" beam-hardening correction.
%|
%| in
%|	fit	strum	initialized by de_ftab_fit()
%|	s1	[N]	s_1 sample values for building table (ftab.sls.sl{1})
%|
%| option
%|	'll'	1	which component (default: 1)
%|	'type'	char	'interp' cubic interpolation (default)
%|			'exp'	log(exp()) form - under development
%|	'show'	0|1	plot it?
%|
%| out
%|	inv1	strum
%|	methods:
%|		inv1.fun(hf1) 	map log values into corrected values
%|			(corresponding to line-integral of material density)
%|		inv1.plot(fit)	show fit
%|
%| Copyright 2008-09-28, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(fit, 'test'), de_ftab_inv1_test, return, end

arg.show = false;
arg.type = 'interp'; % cubic interpolation
arg.ll = 1; % which material, one of {1, ..., LL}
arg = vararg_pair(arg, varargin);

% although we need only the 1D array of s1 values,
% fit.fmfun requires a [n1 1 LL] array.  here the other components are zero.
LL = fit.LL; % usually 1 anyway
sll = s1; sll(1,1,LL) = 0; % [n1 1 L]

% evaluate nonlinear BH function for each of the M spectra
f1 = fit.fmfun(sll); % [n1 1 M]
%f1 = squeeze(f1(:,1,:)); % [n1 M] % changed 2008-12-19
f1 = squeeze(f1); % [n1 M]

st.s1 = s1;
st.f1 = f1;

switch arg.type
case 'interp'
	inv1_fun = @de_ftab_inv1_interp;
case 'exp'
	st = de_ftab_inv1_exp_init(st, fit, arg);
	inv1_fun = @de_ftab_inv1_exp_eval;
otherwise
	fail('unknown type %s', arg.type)
end

meth = {'fun', inv1_fun, '(fh1 [N M]) -> [N M]';
	'plot', @de_ftab_inv1_plot, '(fit)'};

inv1 = strum(st, meth);

if arg.show && im
	inv1.plot(fit);
end


% de_ftab_inv1_exp_init()
function st = de_ftab_inv1_exp_init(st, fit, arg);

s1 = st.s1;
f1 = st.f1;

arg.fit0 = true; % todo
arg.fit0 = false;
arg.chat = true;
arg.thresh = 1e-16;
arg.nexp = 91;

y = exp(s1);
wt = ones(size(y)); % todo
Wh = spdiag(wt, 'nowarn');

MM = size(f1,2);
for mm=1:MM
	f_coef = fit.coef{mm};
	f_exp = fit.mac{mm}(:,1); % for ll=1 only
	emin = 1 / max(f_exp);
	emax = 1 / min(f_exp);
	inv_exp = linspace(emin, emax, arg.nexp)'; % [Ne 1] candidate exps
	Ab = exp(f1(:,mm) * inv_exp'); % [*Ns Ne] "over-complete" basis

	% initial coefficients for each candidate exponent
	if arg.fit0
		fder0 = fit.mac_eff(mm,1);
		inv_der0 = 1 / fder0;
		warg = {'inprodv', inv_exp' / inv_der0};
	else
		warg = {};
	end
	x = wls_simplex(Ab, y(:), Wh, [], warg{:}); % [Ne 1]

	if 1 % reduce to fewer essential terms
		ie = x > arg.thresh; % find key energies
		if arg.chat
			printm('%d coefficients of %d for m=%d', ...
				sum(ie), length(ie), mm)
		end
		inv_exp = inv_exp(ie);

		% final coefficients at key exponents
		if arg.fit0
			fder0 = fit.mac_eff(mm,1);
			inv_der0 = 1 / fder0;
			warg = {'inprodv', inv_exp' / inv_der0}; % shorter now
		else
			warg = {};
		end
		A = Ab(:,ie); % [*Ns Ne] final basis
		x = wls_simplex(A, y(:), Wh, [], warg{:}); % [Ne 1]
	end

	st.exp.coef{mm} = x; % [Ne 1]
	st.exp.exp{mm} = inv_exp;
end


% de_ftab_inv1_exp_eval()
function s1hat = de_ftab_inv1_exp_eval(st, fh1)

MM = size(st.f1,2);
dim = size(fh1);
fh1 = reshapee(fh1, [], MM); % [N M]
s1hat = zeros(size(fh1), 'single');
for mm=1:MM
	tmp = exp(fh1(:,mm) * st.exp.exp{mm}'); % [*N Ne]
	s1hat(:,mm) = log(tmp * st.exp.coef{mm}); % [*N 1]
end
s1hat = reshape(s1hat, dim);

warn 'untested'


% de_ftab_inv1_interp()
% inverse based on interpolation
function s1hat = de_ftab_inv1_interp(st, fh1)

s1 = st.s1;
f1 = st.f1;

MM = size(f1,2);
dim = size(fh1);
fh1 = reshapee(fh1, [], MM);
s1hat = zeros(size(fh1), 'single');
for mm=1:MM
	s1hat(:,mm) = interp1(f1(:,mm), s1, fh1(:,mm), 'pchip', 'extrap');
end
s1hat = reshape(s1hat, dim);


% de_ftab_inv1_plot()
function out = de_ftab_inv1_plot(st, fit)

if nargin ~= 2, fail 'need "fit" argument', end

s1 = st.s1;

s4 = max(s1);
s4 = linspace(-0.0*s4, 1.0*s4, 4*length(s1)+1)'; % [n4] fine sampling
LL = fit.LL;
MM = fit.MM;
sll = s4; sll(1,1,LL) = 0; % [n4 1 L]

f4 = fit.fmfun(sll); % [n4 1 M]
%f4 = squeeze(f4(:,1,:)); % [n4 M]
f4 = squeeze(f4); % [n4 M]

s4i = st.fun(f4); % [n4 M] one correction for each spectrum
err = s4i - repmat(s4, [1 MM]); % [n4 M]

clf, pl = @(n,m) subplot(200 + MM*10 + (n-1)*MM + m);
for mm=1:MM
	pl(1,mm)
	plot(s4, f4(:,mm), 'c-', s4i(:,mm), f4(:,mm), 'y:', ...
		s4, s4 * fit.mac_eff(mm), '--') % monenergetic line
	ir_legend({'true (fine)', 'fit', 'mono'})
	m = num2str(mm);
	axis tight, xlabel 's1', ylabel(['f' m]), title(['inv1 fit m=' m])

	pl(2,mm)
	plot(s4, err(:,mm), '.-')
	axis tight, xlabel 's1', ylabel(['err' m]), title 'inv1 error'
end
prompt

if nargout, out = []; end


% de_ftab_inv1_test()
function de_ftab_inv1_test
stype = 'ps1';
mtype = {'water', 'bone'};
xrs = xray_read_spectra(stype);
mas = xray_read_mac(mtype);
sls = de_ftab_sls;
fm = de_ftab_fm(sls.sll, mas.mac(xrs.en), xrs.Ide);
sl = sls.sl;
fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype, 'show', 0);
itypes = {'exp', 'interp'};
for ii=1:length(itypes)
	itype = itypes{ii};
	inv1 = de_ftab_inv1(fit, sl{1}, 'type', itype);
%	inv1.fun([0 0])
	if im, inv1.plot(fit); end
end
