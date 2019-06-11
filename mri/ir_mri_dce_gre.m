 function sig = ir_mri_dce_gre(r1d, varargin)
%function sig = ir_mri_dce_gre(r1d, varargin)
%|
%| Given dynamic time curves for R1=1/T1, compute corresponding MRI signal
%| for (spoiled) gradient echo (GRE) signal model.
%| It is known this model is inaccurate for imperfect spoiling
%| but it is useful for simulations at least.
%|
%| in
%|	r1d	[Nc Nt]		[1/s] dynamic R1=1/T1 curves
%|				for Nc tissue classes
%|
%| option
%|	'TR'	[1]		[s] repetition time (default: 5e-3)
%|	'TE'	[1]		[s] echo time (default: 0)
%|	'flip'	[1]		[rad] flip angle (default: pi/6)
%|	'M0'	[1] | [Nc]	equilibrium magnetization (default: 1)
%|	'R2s'	[1] | [Nc]	R2* (default: 0)
%|
%| out
%|	sig	[Nc Nt]		[au] signal
%|
%| 2014-08-21 Jeff Fessler and Mai Le, University of Michigan

if nargin == 1 && streq(r1d, 'test'), ir_mri_dce_gre_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.TR = 5e-3; % 5ms TR
arg.TE = 0;
arg.flip = pi/6; % 30 degree flip
arg.M0 = 1;
arg.r2s = 0;
arg = vararg_pair(arg, varargin);

sig = ir_mri_dce_gre_do(r1d, arg.TR, arg.TE, arg.flip, arg.M0, arg.r2s);


% ir_mri_dce_gre_do()
function sig = ir_mri_dce_gre_do(r1d, TR, TE, flip, M0, r2s)

if numel(M0) ~= 1
	if numel(M0) ~= nrow(r1d), fail 'Nc mismatch for M0', end
	M0 = repmat(M0(:), [1 ncol(r1d)]);
end
if numel(r2s) ~= 1
	if numel(r2s) ~= nrow(r1d), fail 'Nc mismatch for R2*', end
end
M = M0 .* exp(-TE * r2s);
if numel(M) > 1
	M = repmat(M(:), [1 ncol(r1d)]);
end
E1 = exp(-TR * r1d);
sig = M .* sin(flip) .* (1 - E1) ./ (1 - E1 .* cos(flip));


% ir_mri_dce_gre_test()
function ir_mri_dce_gre_test

TR = 5e-3; % 5 msec
duration = 4; % [min]
ti = 0:TR:duration;
Ktrans = [3.0 0.6 0.2 0];
kep = [6.0 2.0 1.3 1];
leg = {'rapid', 'moderate', 'slow', 'none'};
Ct = ir_mri_dce_aif1(ti, Ktrans, kep);

f = mri_brainweb_params('grey-matter');
t10 = f.t1 / 1000; % [s]
r10 = 1 / t10; % [1/s]
r10 = [0.9 1.0 1.1 1.0] * r10;

r1d = ir_mri_dce_r1d(r10, Ct);

r2s = 1000 / f.t2s; % R2* in [1/s]
r2s = [0.9 1.0 1.1 1.0]' * r2s;
sig = ir_mri_dce_gre(r1d, 'TR', TR, 'TE', 20e-3, 'r2s', r2s);

if im
	ccoef = corrcoef(sig');
	pr ccoef
	pr ccoef(4)
	pr svd(sig)
	pr cond(sig*sig') % condition number about 2e4

	clf
	plot(ti, sig, '.-')
	legend(leg{:})
	xlabel 't [min]'
	ylabelf 'signal'
end
