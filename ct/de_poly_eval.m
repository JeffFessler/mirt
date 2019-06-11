  function yy = de_poly_eval(xx, coef, expo, varargin)
%|function yy = de_poly_eval(xx, coef, expo, [option])
%|
%| Evaluate multi-dimensional polynomial:
%| y[n,l] = sum_k=1^K coef(k,l) x[n,1]^expo(k,1) * ... * x[n,M]^expo(k,M)
%|
%| in
%|	xx	[(nn) M]	input values
%|	coef	[K L]		coefficients
%|	expo	[K M]		exponents aka powers
%|
%| option
%|	'deriv'	1|0		if 1, then instead return derivative
%|	'basis'	1|0		if 1, then instead return polynomial basis
%|				(coef ignored)
%|
%| out
%|	yy	[(nn) L]	output values
%|		[(nn) L M]	if deriv=1
%|		[(nn) K]	if basis=1
%|
%| Copyright 2008-10-27, Jeff Fessler, University of Michigan

if nargin == 1 && streq(xx, 'test'), de_poly_eval_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.deriv = 0;
arg.basis = false;
arg = vararg_pair(arg, varargin);

dim = size(xx);
MM = dim(end);
dim = dim(1:end-1);

xx = reshapee(xx, [], MM); % [*nn M]

KK = size(expo,1);
if size(expo,2) ~= MM
	fail 'bad expo size'
end

if arg.basis
	if ~isempty(coef) || arg.deriv
		fail 'bad input: need empty coef and deriv=0'
	end
	yy = de_poly_eval_basis(xx, expo);
	yy = reshape(yy, [dim KK]);
return
end

LL = size(coef,2);
if KK ~= size(coef,1)
	fail 'bad coef size'
end

switch arg.deriv
case 0
	yy = de_poly_eval_work(xx, coef, expo);
	yy = reshape(yy, [dim LL]);

case 1
	yy = zeros(prod(dim), LL, MM); % [*nn L M]
	for mm=1:MM
		expom = expo; % [K M]
		coefm = coef; % [K L]
		expom(:,mm) = expo(:,mm) - 1;

		for kk=1:KK
			if expo(kk,mm) == 0
				expom(kk,mm) = 0;
				coefm(kk,:) = 0;
			else
				coef(kk,:) = coef(kk,:) .* expo(kk,mm);
			end
		end
		yy(:,:,mm) = de_poly_eval_work(xx, coefm, expom);
	end
	yy = reshape(yy, [dim LL MM]);

otherwise
	fail 'only deriv = 0 or 1 done'
end


%
% de_poly_eval_work()
%
function yy = de_poly_eval_work(xx, coef, expo)

KK = size(expo,1);
LL = size(coef,2);
nn = size(xx,1);
yy = zeros(nn, LL); % [*nn L]
for kk=1:KK
	po = expo(kk,:); % [1 M] powers
	co = coef(kk,:); % [1 L] coefficients

	po = repmat(po, [nn 1]); % [*nn M]
	tmp = xx .^ po; % [*nn M]
	tmp = prod(tmp, 2); % [*nn 1]
	yy = yy + tmp * co;
end


%
%function basis = de_poly_eval_basis(xx, expo)
%
% Basis for multi-dimensional polynomial:
% basis[n,k] = x[n,1]^expo(k,1) ... x[n,M]^expo(k,M)
%
% in
%	xx	[N M]		input values
%	expo	[K M]		exponents aka powers
%
% out
%	basis	[N K]		basis vectors
%
function basis = de_poly_eval_basis(xx, expo)

[NN MM] = size(xx);
[KK MM] = size(expo);

basis = zeros(NN, KK); % [N K]
for kk=1:KK
	po = expo(kk,:); % [1 M] powers
	po = repmat(po, [NN 1]); % [N M]
	tmp = xx .^ po; % [N M]
	basis(:,kk) = prod(tmp, 2); % [N 1]
end


%
% de_poly_eval_test()
%
function de_poly_eval_test

xx = ndgrid_jf('mat', 0:4, 0:5);
expo = [0 0; 1 0; 0 1; 1 1];
coef = [1 2 3; 0 0 1; 1 1 0; 2 5 6];
yy = de_poly_eval(xx, coef, expo);

x1 = xx(:,:,1);
x2 = xx(:,:,2);
ll = size(coef,2);
co = coef(:,ll);
tmp = co(1) + co(2) * x1 + co(3) * x2 + co(4) * x1 .* x2;
jf_equal(yy(:,:,ll), tmp)


yy = de_poly_eval(xx, coef, expo, 'deriv', 1);
tmp = 0*co(1) + co(2) + 0 * co(3) + co(4) * x2;
jf_equal(yy(:,:,3,1), tmp)

tmp = 0*co(1) + 0 * co(2) + co(3) + co(4) * x1;
jf_equal(yy(:,:,3,2), tmp)

yy = de_poly_eval(xx, [], expo, 'basis', 1);
im(yy)
