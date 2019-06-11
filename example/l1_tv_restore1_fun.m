 function xs = l1_tv_restore1_fun(yi, A, varargin)
%function xs = l1_tv_restore1_fun(yi, A, varargin)
%|
%| Robust image restoration using a l1 data fit term and a TV-like regularizer.
%|
%| see l1_tv_restore1.m
%| min_x |yi - A(ti) * x|_1 + beta (|Cx x|_1 + |Cy x|_1)
%| where |.|_1 is approximated by a hyperbola
%|
%|
%| in
%|	yi	data
%|	A	system matrix
%|
%| option
%|	several - see code ...
%|
%| out
%|	xs	estimate(s)
%|
%| Copyright 2010-05-24, Jeff Fessler, University of Michigan

if nargin == 1 && streq(yi, 'test')
	l1_tv_restore1
	clear xs
return
end

if nargin < 2, help(mfilename), error(mfilename), end

arg.niter = 20;
arg.delta1 = 0.2; % round corner of |t| in data fit to approximate by hyperbola
arg.delta2 = 0.2; % round corner of |t| in penalty
arg.xinit = [];
arg.mask = [];
arg.l2b = -6;
arg.curvtype = 'oc';
arg.args = {};
arg = vararg_pair(arg, varargin);

if isempty(arg.xinit)
	arg.xinit = yi; % assume restoration
end

if isempty(arg.mask)
	arg.mask = true(size(yi));
end

% R = R_null; % null regularizer 
R = Reg1(arg.mask, 'type_denom', 'matlab', ...
	'beta', 2^arg.l2b, 'pot_arg', {'hyper2', arg.delta2});

data = {yi(:), arg.delta1};
xs = pl_pcg_qs_ls(arg.xinit(:), A, data, ...
	@l1_tv_dercurv, R, 'niter', arg.niter, ...
	'curvtype', arg.curvtype, ...
	arg.args{:});


% l1_tv_dercurv()
function [deriv curv] = l1_tv_dercurv(data, yp, curvtype)
yi = data{1};
delta = data{2};
pot = potential_fun('hyper2', delta);
t = yp - yi;
switch curvtype
case 'pc'
	curv = 1;
	deriv = pot.dpot(t);
case 'oc'
	curv = pot.wpot(t);
	deriv = curv .* t;
otherwise
	fail('curvtype %s', curvtype)
end


% R_null()
function R = R_null;
R.dercurv = @R_null_dercurv;
R.C1 = 0;

% R_null_dercurv()
function [a b] = R_null_dercurv(arg1, arg2);
a = 0;
b = 0;
