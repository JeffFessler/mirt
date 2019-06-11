 function r1d = ir_mri_dce_r1d(r10, Ct, varargin)
%function r1d = ir_mri_dce_r1d(r10, Ct, varargin)
%|
%| Make time series of R1=1/T1 values given baseline value(s) R10
%| and contrast agent time-concentration curves Ct.
%|
%| in
%|	r10	[1] | [Nc]	[1/s] baseline R1 values
%|				for Nc tissue classes
%|	Ct	[Nc Nt]		[mMmol] tissue constrast concentration curves
%|				for Nt time points
%|
%| option
%|	'r1'	[1]		relaxivity [1/mMol 1/s]; default: 4.5 (for 3T)
%|
%| out
%|	r1d	[Nc Nt]		[1/s] dynamic R1=1/T1 curves
%|
%| 2014-08-21 Jeff Fessler and Mai Le, University of Michigan

if nargin == 1 && streq(r10, 'test'), ir_mri_dce_r1d_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.r1 = 4.5; % [mMol^-1 sec^-1] for Gd at 3T, Table 2 of sasaki:05:eea
arg = vararg_pair(arg, varargin);

if numel(r10) ~= 1 && numel(r10) ~= nrow(Ct)
	fail 'Nc mismatch'
end
r1d = ir_mri_dce_r1d_do(r10(:), Ct, arg.r1);


function r1d = ir_mri_dce_r1d_do(r10, Ct, r1)

if numel(r10) ~= 1
	r10 = repmat(r10, [1 ncol(Ct)]);
end
r1d = r10 + r1 * Ct;


% ir_mri_dce_r1d_test()
function ir_mri_dce_r1d_test

TR = 5e-3; % 5 msec
duration = 4; % [min]
ti = linspace(0, duration, 401); % [min]
Ktrans = [3.0 0.6 0.2 0];
kep = [6.0 2.0 1.3 1];
leg = {'rapid', 'moderate', 'slow', 'none'};
Ct = ir_mri_dce_aif1(ti, Ktrans, kep);

f = mri_brainweb_params('grey-matter');
t10 = f.t1 / 1000; % [s]
r10 = 1 / t10; % [1/s]

r1d = ir_mri_dce_r1d(r10, Ct);
if im
	clf
	subplot(121)
	plot(ti, r1d, '.-'), ylabelf 'R1 [1/s]'
	legend(leg{:})
	t1d = 1 ./ r1d;
	subplot(122)
	plot(ti, 1000*t1d, '.-'), ylabelf 'T1 [ms]'
	xlabel 't [min]'
end
