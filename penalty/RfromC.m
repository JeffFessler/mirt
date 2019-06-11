 function R = RfromC(varargin)
%function R = RfromC(varargin)
%|
%| Create regularization object from "analysis" matrix (or fatrix2) C1.
%|
%| R(x) = sum_k wt_k pot([C1*x]_k)
%|
%| option
%|	'C1'		from Godwt1 or such.  default: 1
%|			must be capable of providing abs(C1)
%|	'pot'		from potential_fun. default: quadratic
%|	'wt'		size of C1*x, including beta. default: 1
%|
%| out
%|	R		strum with methods:
%|			R.cgrad(R,x) column gradient
%|			R.penal(R,x) penalty value
%|			R.denom(R,x) denominator for SQS
%|			R.C(R,x) - of limited use.  C = diag(wt) * C1
%|
%| Copyright 2011-1-21, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(varargin{1}, 'test'), RfromC_test, return, end

arg.C1 = 1;
arg.pot = potential_fun('quad');
arg.wt = 1;
arg = vararg_pair(arg, varargin);

meth = {
	'cgrad', @RfromC_cgrad, '(R,x)';
	'penal', @RfromC_penal, '(R,x)';
	'denom', @RfromC_denom, '(R,x)';
	'C', @RfromC_C, '()';
};

R = strum(arg, meth);


% RfromC_penal()
function penal = RfromC_penal(dummy, arg, x)
tmp = arg.C1 * x;
pot = arg.pot.potk(tmp);
penal = sum(arg.wt(:) .* pot(:));


% RfromC_cgrad()
function cgrad = RfromC_cgrad(dummy, arg, x)
tmp = arg.C1 * x;
wpot = arg.pot.wpot(tmp);
if numel(arg.wt) == 1
	wt = wpot * arg.wt; % scalar wt
else
	wt = wpot .* reshape(arg.wt, size(wpot));
end
cgrad = arg.C1' * (wt .* tmp);


% RfromC_denom()
function denom = RfromC_denom(dummy, arg, x)
Ca = abs(arg.C1);
tmp = Ca * x;
wpot = arg.pot.wpot(tmp);
if numel(arg.wt) == 1
	wt = wpot * arg.wt; % scalar wt
else
	wt = wpot .* reshape(arg.wt, size(wpot));
end
c1 = Ca * ones(size(x));
denom = Ca' * (wt .* c1);


% RfromC_C()
function C = RfromC_C(dummy, R)
D = Gdiag(sqrt(arg.wt), 'mask', true(size(arg.wt)));
C = D * arg.C1;


function RfromC_test
mask = true(2^5, 2^4); mask(2) = false; % stress
Ch = Gblur(mask, 'psf', [0; 1; -1], 'type', 'fft,same');
Cv = Gblur(mask, 'psf', [0 1 -1], 'type', 'fft,same');
C1 = [Ch; Cv];
%C1 = Godwt1(mask); % todo: need "abs"
pot = potential_fun('lange1', 0.01);
wt = 3;
R = RfromC('pot', pot, 'C1', C1, 'wt', 3);
x = zeros(size(mask)); x(end/2+1,end/2+1)=1;
cx = R.C1 * x;
im plc 2 2
cgrad = R.cgrad(R, x);
denom = R.denom(R, x);
penal = R.penal(R, x);
im(1, x), cbar
im(2, cx), cbar
im(3, cgrad), cbar
im(4, denom), cbar
