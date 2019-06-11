 function f = mri_brainweb_params(label, varargin)
%function f = mri_brainweb_params(label, varargin)
%|
%| tissue parameters for brainweb brain images, based on
%| http://mouldy.bic.mni.mcgill.ca/brainweb/tissue_mr_parameters.txt
%| units of t* are msec
%|
%| in
%|	label		{0,1,...,10} or 'fat', etc.
%| option
%|	b0		1.5 (default) or 3 (todo) (in Tesla)
%| out
%|	f.t1,t2,t2s,pd,label,name
%|
%| 2012-05-15, Jeff Fessler, Univ. of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(label, 'test'), mri_brainweb_params_test, return, end

arg.b0 = 1.5;
arg = vararg_pair(arg, varargin);

if arg.b0 ~= 1.5
	fail('only 1.5T done now')
end

if ischar(label)
	for ii=0:10
		f = mri_brainweb_params(ii);
		if streq(f.name, label)
			return
		end
	end
	fail('label "%s" not found', label)
return
end

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

f.label = label;


% mri_brainweb_params_test()
function mri_brainweb_params_test
for ii=0:10
	f = mri_brainweb_params(ii);
	f.name; f.t1; f.t2; f.t2s; f.pd;
end
