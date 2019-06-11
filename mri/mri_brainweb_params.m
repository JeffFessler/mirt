 function f = mri_brainweb_params(label, varargin)
%function f = mri_brainweb_params(label, varargin)
%|
%| Tissue parameters for BrainWeb synthetic brain images:
%| 1.5T values from BrainWeb database:
%|	http://mouldy.bic.mni.mcgill.ca/brainweb/tissue_mr_parameters.txt
%| 3.0T T1, T2 values taken from literature when possible (citations below)
%|	[1] Wansapura et al. "NMR Relaxation Times in the Human Brain at
%|		3.0 Tesla." JMRI vol. 9, no. 4, Apr. 1999.
%|	[2] Stanisz et al. "T1, T2 Relaxation and Magnetization Transfer in
%|		Tissue at 3T." MRM vol. 54, no. 3, Sept. 2005.
%|
%|	Otherwise, T1, T2 estimated as twice the 1.5T value.
%|	M0 values copied from 1.5T value
%|	T2* estimated as if R2prime_3T = 2*R2prime_1.5T
%|
%| Units of t* are msec
%|
%| in
%|	label		{0,1,...,10} or 'fat', etc.
%| option
%|	b0		1.5 (default) or 3 (Tesla)
%| out
%|	f.t1,t2,t2s,pd,label,name
%|
%| 1.5T version:	2012-05-15, Jeff Fessler, Univ. of Michigan
%| 3.0T extension:	2015-05-11, Gopal Nataraj, Univ. of Michigan

if nargin < 1, ir_usage, end
if streq(label, 'test'), mri_brainweb_params_test, return, end

arg.b0 = 1.5; % default
arg = vararg_pair(arg, varargin);

if arg.b0 ~= 1.5
	warn('3T values partially extrapolated')
end

% If the input is a string, run with the numerical equivalent
if ischar(label)
	for ii = 0:10
		f = mri_brainweb_params(ii, varargin{:});
		if streq(f.name, label)
			return
		end
	end
	fail('label "%s" not found', label)
return
end

switch arg.b0
case 1.5
	f = mri_brainweb_params_1p5(label);
case 3
	f = mri_brainweb_params_3p0(label);
otherwise
	fail('Values for field strength %g unknown', arg.b0)
end

end % mri_brainweb_params()


% f = mri_brainweb_params_1p5()
function f = mri_brainweb_params_1p5(label)


switch label
case 0
	f.name = 'background';
	f.t1 = 0;
	f.t2 = 0;
	f.t2s = 0;
	f.pd = 0;

case 1
	f.name = 'csf';
	f.t1 = 2569;
	f.t2 = 329;
	f.t2s = 58;
	f.pd = 1;

case 2
	f.name = 'grey-matter';
	f.t1 = 833;
	f.t2 = 83;
	f.t2s = 69;
	f.pd = 0.86;

case 3
	f.name = 'white-matter';
	f.t1 = 500;
	f.t2 = 70;
	f.t2s = 61;
	f.pd = 0.77;

case 4
	f.name = 'fat';
	f.t1 = 350;
	f.t2 = 70;
	f.t2s = 58;
	f.pd = 1;

case 5
	f.name = 'muscle/skin';
	f.t1 = 900;
	f.t2 = 47;
	f.t2s = 30;
	f.pd = 1.0;

case 6
	f.name = 'skin';
	f.t1 = 2569;
	f.t2 = 329;
	f.t2s = 58;
	f.pd = 1;

case 7
	f.name = 'skull';
	f.t1 = 0;
	f.t2 = 0;
	f.t2s = 0;
	f.pd = 0;

case 8
	f.name = 'glial matter';
	f.t1 = 833;
	f.t2 = 83;
	f.t2s = 69;
	f.pd = 0.86;

case 9
	f.name = 'meat';
	f.t1 = 500;
	f.t2 = 70;
	f.t2s = 61;
	f.pd = 0.77;

case 10
	f.name = 'ms-lesion';
	f.t1 = 752;
	f.t2 = 237;
	f.t2s = 204;
	f.pd = 0.76;

otherwise
	fail('bad label %d', label)
end
end % mri_brainweb_params_1p5()


% mri_brainweb_params_3p0()
function f = mri_brainweb_params_3p0(label)

% Several values estimated from same label, but at 1.5T
f_1p5T = mri_brainweb_params_1p5(label);
f.name = f_1p5T.name; % Copy name over from 1.5T

switch label
case 0 % Background: estimated from 1.5T
	f.t1 = t1_3T_est(f_1p5T);
	f.t2 = t2_3T_est(f_1p5T);

case 1 % CSF: not much info, estimate order of magnitude
	f.t1 = 3000;
	f.t2 = 350;

case 2 % Grey matter: from [1]
	f.t1 = 1331;
	f.t2 = 110;

case 3 % White matter: from [1]
	f.t1 = 832;
	f.t2 = 79.6;

case 4 % Fat: estimated from 1.5T
	f.t1 = t1_3T_est(f_1p5T);
	f.t2 = t2_3T_est(f_1p5T);

case 5 % Muscle/skin: from [2]
	f.t1 = 1412;
	f.t2 = 50;

case 6 % Skin: same as CSF
	f.t1 = 3000;
	f.t2 = 350;

case 7 % Skull: estimated from 1.5T
	f.t1 = t1_3T_est(f_1p5T);
	f.t2 = t2_3T_est(f_1p5T);

case 8 % Glial matter: same as grey matter
	f.t1 = 1331;
	f.t2 = 110;

case 9 % Meat: same as white matter
	f.t1 = 832;
	f.t2 = 79.6;

% MS-Lesion: unable to find well-cited, agreed-upon values
% T2* also likely to change differently due to iron overload
% Safe bet: extrapolated from 1.5T values.
case 10
	f.t1 = t1_3T_est(f_1p5T);
	f.t2 = t2_3T_est(f_1p5T);

otherwise
	fail('bad label %d', label)
end

% t2s, pd always estimated from 1.5T
f.t2s = t2s_3T_est(f_1p5T, f.t2);
f.pd = pd_3T_est(f_1p5T);
f.label = label;
end % mri_brainweb_params_3p0()

% Estimating t1 at 3T from t1 at 1.5T
function t1_3T = t1_3T_est(f_1p5T)
	relax_factor = sqrt(3/1.5); % sqrt relation (?)
	t1_3T = relax_factor * f_1p5T.t1;
end

% Estimating t2 at 3T from t2 at 1.5T
function t2_3T = t2_3T_est(f_1p5T)
	relax_factor = sqrt(3/1.5); % sqrt relation (?)
	t2_3T = relax_factor * f_1p5T.t2;
end


% Estimating t2s at 3T from t2s at 1.5T
function t2s_3T = t2s_3T_est(f_1p5T, t2_3T)
	R2p_1p5T = div0(1, f_1p5T.t2s) - div0(1, f_1p5T.t2); % R2p = R2s - R2
	R2p_3T = 2 * R2p_1p5T; % R2p_3T doubled
	R2s_3T = div0(1, t2_3T) + R2p_3T;
	t2s_3T = div0(1, R2s_3T);
end

% Estimating pd at 3T from pd at 1.5T
function pd_3T = pd_3T_est(f_1p5T)
	pd_3T = f_1p5T.pd; % Just equated for now
end


% mri_brainweb_params_test()
function mri_brainweb_params_test
for b0 = [1.5 3]
	for ii = 0:10
		f = mri_brainweb_params(ii, 'b0', b0);
		if im
			f.name; f.t1; f.t2; f.t2s; f.pd;
			disp(f);
		end
	end
end
end % mri_brainweb_params_test()
