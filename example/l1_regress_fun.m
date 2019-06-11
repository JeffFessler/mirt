 function [x1, x2] = l1_regress_fun(ti, yi, varargin)
%function [x1, x2] = l1_regress_fun(ti, yi, varargin)
%|
%| see l1_regress_example.m
%| min_x |yi - A(ti) * x|_1
%| where |.|_1 is approximated by a hyperbola
%|
%| Copyright 2005-4-24, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ti, 'test')
	l1_regress_example
	clear x1 x2
return
end
if nargin < 2, help(mfilename), error(mfilename), end

arg.niter = 100;
arg.delta = 0.2; % round corner of |t| to approximate by hyperbola
arg.linear = false;
arg.A = [];
arg.xinit = [];
arg.isave = 'all';
%arg.l2b = -6;
arg = vararg_pair(arg, varargin);

if isempty(arg.A)
	if arg.linear
		arg.A = [ti];
	else
		arg.A = [ones(size(ti)) ti];
	end
%	x2 = regress(yi, arg.A); % ordinary l_2 regression
	x2 = arg.A \ yi; % ordinary l_2 regression
	if isempty(arg.xinit)
		arg.xinit = x2;
	end
end

if 1 % new way
	% null regularizer 
	R = R_null;
	x1 = pl_pcg_qs_ls(arg.xinit, arg.A, {yi, arg.delta}, @l1_dercurv, ...
		R_null, 'niter', arg.niter, 'isave', arg.isave);
return
end

%
% below here is old way - do not use!
%

%{
G = sparse(zeros(1, 2));
G = Gsparse(G, 'idim', [2 1], 'odim', [1 1]);

% trick: this roughness penalty will serve as the data-fit term! 
R = Rreg1([1; 1], 'type_denom', 'matlab', ...
	'offsets', 0, ... % trick for identity
	'beta', 2^(-arg.l2b), ... % trick: negative because likelihood! 
	'pot_arg', {'hyper3', arg.delta});

R.pot = potential_shift(R.pot, yi); % data!
R.C1 = [ones(size(ti)) ti]; % predictor matrix, aka, system matrix!
R.E = ( diag(sum(abs(R.C1'))) * abs(R.C1) )';
R.denom = @(R, x) R.E * (R.C1 * x);

xinit = x2;
x1 = pwls_sps_os(x2, 0, [], G, R, arg.niter, [-inf inf], [], [], 1);
%}


% l1_dercurv()
function [deriv, curv] = l1_dercurv(data, yp, curvtype)
yi = data{1};
delta = data{2};
jf_equal(curvtype, 'pc')
%curv = 1;
pot = potential_fun('hyper2', delta);
t = yp - yi;
curv = pot.wpot(t);
deriv = curv .* t;


% R_null()
function R = R_null
R.dercurv = @R_null_dercurv;
R.C1 = 0;

% R_null_dercurv()
function [a, b] = R_null_dercurv(arg1, arg2)
a = 0;
b = 0;
