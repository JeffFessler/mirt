  function sh = de_ftab_s_iter(fit, fh, varargin)
%|function sh = de_ftab_s_iter(fit, fh, [options])
%| estimate s from fh by iterative LS
%| sh = argmin{s} | fh - fm(s) |^2
%| in
%|	fit			from de_ftab_fit() / de_ftab_curv()
%|	fh	[(Nd) M]	estimates of f (nonlinear BH function)
%| option
%|	'niter'			# of iterations (default: 1 - not enough!)
%|				(if 0, then just initial estimate returned)
%|	'init'	[(Nd) L]	initial estimates (default: linear inverse)
%|	'ctype'	char		curvature type (see de_ftab_curv.m)
%|				default '': inherit from de_ftab_curv.m
%| out
%|	sh	[(Nd) L]	estimates of s (component density integrals)
%|
%| Copyright 2006-05-21, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fit, 'test'), de_ftab_s_iter_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.ctype = '';
arg.niter = 1;
arg.init = [];
arg.strue = []; % optional for testing
arg.nprint = 100; % every 100th iteration show rms
arg = vararg_pair(arg, varargin);

LL = fit.LL;
Nd = size(fh); MM = Nd(end); Nd = Nd(1:end-1);
fh = reshapee(fh, [], MM); % [*Nd M]
if MM ~= fit.MM, error 'M', end

if ~isempty(arg.strue)
	arg.strue = reshapee(arg.strue, prod(Nd), LL); % [*Nd L]
end

% initialize
if isempty(arg.init)
	sh = (pinv(fit.mac_eff) * fh')'; % linear inverse (perfect if mono)

	if ~isempty(arg.strue)
		pr 'rms(sh - arg.strue)' % report initial error, for testing!
	end
else
	sh = max(arg.init, 0);
	sh = reshape(sh, prod(Nd), LL); % [*Nd L]
end

if arg.niter < 1, return, end

% curvature
switch arg.ctype
case ''
%	fstep = @(sh, fm, fh) 1 ./ sum(fit.curv2(sh, fm, fh), 2); % [*Nd 1]
	fstep = @(sh, fm, fh) 1 ./ sum(fit.ls_curv(sh, fh, fm), 2); % [*Nd 1]

case 'old'
	[curv1 curv2] = de_ftab_curv_ub(fit); % each is [MM 1]
	fcurv2 = @(fm, fh) max(fh - fm, 0) * curv2; % [*Nd M] * [M 1] -> [*Nd 1]
%	fcurv2 = @(fm, fh) max(fh, 0) * curv2; % [*Nd M] * [M 1] -> [*Nd 1] TOO SLOW!?
	fstep = @(sh, fm, fh) 1 ./ (sum(curv1) + fcurv2(fm, fh)); % [*Nd 1]
	step = fstep(0, 0, fh); % [*Nd 1] precomputed!
	minmax(step)
	step = repmat(step, [1 LL]);

otherwise
	fstep = @(sh, fm, fh) 1 ./ ...
		sum(fit.ls_curv(sh, fh, fm, 'ctype', arg.ctype), 2); % [*Nd 1]
%	fail('bad ctype: %s', arg.ctype)
end

ticker reset
for ii=1:arg.niter % iterate
	ticker(mfilename, ii, arg.niter)

	fm = fit.fmfun(sh); % [*Nd M]
%	cost = mean(col(fh - fm).^2);

	if 1
		fgrad = fit.fgrad(sh); % [*Nd L M]
		tmp = repmat(fh - fm, [1 1 LL]); % [*Nd M L]
		tmp = permute(tmp, [1 3 2]); % [*Nd L M]
		fgrad = fgrad .* tmp;
	else
%		fgrad = fit.ls_grad(fh, sh);
	end

	step = fstep(sh, fm, fh); % [*Nd 1]
	step = repmat(step, [1 LL]); % [*Nd L]
%	minmax(step)
	sh = sh + step .* sum(fgrad, 3);

	if ~isempty(arg.strue) && (~rem(ii, arg.nprint) || ii==arg.niter)
		pr 'rms(sh - arg.strue)' % report error, for testing!
	end
end

sh = reshapee(sh, Nd, LL); % [(Nd) L]


% de_ftab_curv_ub()
% curv* are [M 1] upper bounds
%
function [curv1, curv2] = de_ftab_curv_ub(fit)

MM = fit.MM;
new = fit.fgrad(zeros([1 fit.LL])); % [1 L M]
curv1 = squeeze(sum(new.^2,2)); % [1 1 M] -> [M 1]

neghess = -fit.fhess(zeros([1 fit.LL])); % [1 L L M]
curv2 = zeros(MM,1);
for mm=1:MM
	h0 = squeeze(neghess(1,:,:,mm)); % [1 L L 1] -> [L L]
	curv2(mm) = norm(h0);
end

return

% old way:

oldcurv1 = zeros(MM,1);
oldcurv2 = zeros(MM,1);

for mm=1:MM
	alf = fit.coef{mm}; % [ne 1]
	mac = fit.mac{mm}; % [ne L]
	g0 = mac' * alf; % gradient of fit at s=0 (largest point) 
	h0 = mac' * diag(alf) * mac - g0 * g0'; % hessian at s=0 (largest?)
	oldcurv1(mm) = norm(g0)^2;
	oldcurv2(mm) = norm(h0);
end

equivs(curv1, oldcurv1)
equivs(curv2, oldcurv2)


% de_ftab_s_iter_test1
% L=1 M=1 case for sanity checking
%
function de_ftab_s_iter_test1
stype = 'ps1';
stype = 'poly1,100';
xrs = xray_read_spectra(stype);
mtype = 'water';
mas = xray_read_mac(mtype);
s1 = linspace(0, 50, 26)';
fm = de_ftab_fm(s1, mas.mac(xrs.en), xrs.Ide);
fit = de_ftab_fit({s1}, fm, 'type', 'exp', 'mtype', mtype);
ctype = 'newt';
%ctype = 'pre10'; % slow!
fit = de_ftab_curv(fit, 'ctype', ctype);

%fm = fit.fmfun(s1); % test with model fm, not original

% iterate
sh = de_ftab_s_iter(fit, fm, 'strue', s1, 'ctype', '', ...
	'niter', 5, 'nprint', 1);

if im
	clf, plot(s1, s1, 'c-', s1, sh, 'y.-')
	grid, axis equal, axis square
	xlabel 'true', ylabel 'noiseless estimate'
	title(fit.ctype)
prompt
end


% de_ftab_s_iter_test2
%
function de_ftab_s_iter_test2
sl{1} = linspace(0, 50, 26);
sl{2} = linspace(0, 30, 31);
xrs = xray_read_spectra('ps1');
mtype = {'water', 'bone'};
mtype = {'water', 'aluminum'};
mas = xray_read_mac(mtype);
sll = ndgrid_jf('mat', sl{:});
fm = de_ftab_fm(sll, mas.mac(xrs.en), xrs.Ide);
fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype);
fit = de_ftab_curv(fit, 'ctype', 'newt');
%fit = de_ftab_curv(fit, 'ctype', 'pre2');

fm = fit.fmfun(sll); % test with model fm, not original

if 1 % initialize with polynomial inverse
	mac = xray_make_mac(xrs, mas);
	T = pinv(mac.bar);
	inv2 = de_ftab_inv2(fit, sl, 'T', T);
	s0 = inv2.fun(fm);
	s0 = reshapee(s0, [], 2);
else
	% initializing with linear inv.
	s0 = de_ftab_s_iter(fit, fm, 'niter', 0);
end

strue = reshapee(sll, [], 2); % true
pr rms(s0 - strue) % initial error - large!

if 1 % picture of error of initial guess vs (s_1,s_2)
	pr fit.mac_eff
	pr cond(fit.mac_eff)
	tmp = reshape(sqrt(mean((s0 - strue).^2, 2)), size(sll(:,:,1)));
	if im
		clf, im(sl{:}, tmp), cbar % error is smallest at (0,0)
		xlabel(mtype{1}), ylabel(mtype{2})
	end
prompt
end

% todo: run profiler, or pre-tabulate the inverse...
% iterate; initializing with linear inv.
sh = de_ftab_s_iter(fit, fm, 'init', s0, 'niter', 90, ...
		'strue', sll, 'nprint', 10);
pr 'rms(reshapee(sh, [], 2) - strue)'

%sh = reshape(sh, size(sll));
if im
	clf, plot(sl{1}, sh(:,:,1), 'c', sl{2}, sh(:,:,2), 'y')
	grid, axis equal, axis square
	xlabel 'true', ylabel 'noiseless estimate'
prompt
%	plot(fm(:,:,1), fm(:,:,2), '.')
end


% de_ftab_s_iter_test
%
function de_ftab_s_iter_test
de_ftab_s_iter_test2
de_ftab_s_iter_test1
