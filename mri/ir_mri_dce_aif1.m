 function [Ct, Cp, Ca, ht] = ir_mri_dce_aif1(ti, Ktrans, kep, varargin)
%function [Ct, Cp, Ca, ht] = ir_mri_dce_aif1(ti, Ktrans, kep, varargin)
%|
%| Generate tissue contrast concentration curves according to Yang et al.
%| "Multiple reference tissue..." MRM Dec. 07 (yang:07:mrt)
%|
%| in
%|	ti		[Nt]		[min] time samples
%|	Ktrans		[Nc]		[1/min] for Nc tissue classes
%|	kep		[Nc]		[1/min] ""
%|
%| option
%|	'dt'		[1]		[min] inner fine time spacing
%|						(default: 1/60 min)
%|	't0'		[1] | [Nc]	[min] injection wait + transport time(s)
%|						(default: 0.25 min = 15s)
%|	'h_power'	[1]		(default: 4) after [11] in yang:07:mrt
%|	'h_beta'	[1]		(default: 0.03 [min]) in Fig. 2 of ""
%|	'chat'		0|1		verbose? (default: 0)
%|
%| out
%| 	Ct		[Nc Nt]		[mMol] tissue concentration curves
%|	Cp		[Nc Nt]		[mMol] local plasma input conc. (AIF)
%| 	Ca 		[1 Nt]		[mMol] aorta (rarely needed)
%| 	ht 		[Nc Nt]		[?] dispersion function (rarely needed)
%|
%| For more, see p.82 of Khalsa thesis.
%|
%| 2014-03-19 Mai Le
%| 2014-08-21 Jeff Fessler

if nargin == 1 && streq(ti, 'test'), ir_mri_dce_aif1_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.h_power = 4; % [unitless] after [11] in yang:07:mrt
arg.h_beta = 0.03; % [min] in Fig. 2 of yang:07:mrt
arg.t0 = 0.25; % injection wait + bolus arrival time delay [min]
arg.dt = 1/60; % [minutes] fine sampling for convolution
arg.chat = 0;
arg = vararg_pair(arg, varargin);

assert(all(Ktrans(:) >= 0), 'invalid Ktrans value')
assert(all(kep(:) >= 0), 'invalid kep value')

Nt = numel(ti);
Nc = numel(Ktrans);
jf_equal(size(Ktrans), size(kep))

Ct = zeros(Nc, Nt, 'single');
Cp = zeros(Nc, Nt, 'single');
ht = zeros(Nc, Nt, 'single');

t0 = arg.t0;
if numel(t0) > 1
	jf_equal(size(arg.t0), size(Ktrans))
else
	t0 = repmat(t0, [1 Nc]);
end

for ic=1:Nc
	[Ct(ic,:) Cp(ic,:) Ca ht(ic,:)] = ir_mri_dce_aif1_do(ti(:)', ...
		Ktrans(ic), kep(ic), ...
		arg.dt, arg.h_power, arg.h_beta, t0(ic), arg.chat);
end


% ir_mri_dce_aif1_do()
% time curves for a single choice of Ktrans and kep
function [Ct Cp Ca ht] = ir_mri_dce_aif1_do(ti, Ktrans, kep, ...
				dt, h_power, h_beta, t0, chat)

t = 0:dt:max(ti(:)); % [minutes] fine sampled

% aorta AIF based on eqn. [12] of yang:07:mrt
Ca = 7.5527 * exp(-(t-0.171).^2 / 0.00605) ...
	+ 1.0003 * exp(-(t-0.364).^2 / 0.035912) ...
	+ 1.064 * exp(-0.083*t) ./ (1 + exp(-37.772*(t-0.482)));

% construct transport function h(t) based on eqn. [11] of yang:07:mrt
% that eqn seems to be missing the "t >= t0" condition but we include it
b = h_beta; % minutes
a = h_power; % unitless power
ht = (b^(-a) .* ((t-t0).^(a-1)) .* exp(-(t-t0)/b) / gamma(a)) .* (t >= t0);

% local plasma input function based on eqn. [10] of yang:07:mrt
Cp = dt * ir_conv1_trim(Ca, ht);

% tissue concentration for "simple Tofts model" based on eqn. [2] of yang:07:mrt
Ct = dt * ir_conv1_trim(Ktrans * Cp, exp(-kep*t));

% may not use full duration, throw away leftover
Ca = interp1(t, Ca, ti, 'linear', 'extrap');
Cp = interp1(t, Cp, ti, 'linear', 'extrap');
Ct = interp1(t, Ct, ti, 'linear', 'extrap');
ht = interp1(t, ht, ti, 'linear', 'extrap');

assert(all(Ca(:) >= 0),'Ca should not have negative values!')
assert(all(Cp(:) >= 0),'Cp should not have negative values!')
assert(all(Ct(:) >= 0),'Ct should not have negative values!')


% ir_conv1_trim()
% akin to 'same' option but trims from tail end rather than start
function y = ir_conv1_trim(x, h)
y = conv(x, h);
y = y(1:length(x));


% ir_mri_dce_aif1_test
function ir_mri_dce_aif1_test

duration = 4; % [min]
ti = linspace(0, duration, 401);
Ktrans = [3.0 0.6 0.2 0];
kep = [6.0 2.0 1.3 0];
leg = {'rapid', 'moderate', 'slow', 'none'};
[Ct Cp Ca h] = ir_mri_dce_aif1(ti, Ktrans, kep, 't0', [0.5 0.6 0.7 0]);

if im
	clf
	subplot(221), plot(ti, Ca, '.-'), ylabel 'C_a(t) aorta'
	subplot(222), plot(ti, h, '.-'), ylabel 'h(t)'
	subplot(223), plot(ti, Cp, '.-'), ylabel 'C_p(t) plasma local AIF'
	subplot(224), plot(ti, Ct, '.-'), ylabel 'C_t(t) tissue curve [mmol]'
	legend(leg{:})
	xlabelf 't [min]'
end
