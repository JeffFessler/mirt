  function ob = diag_sp(diag, varargin)
%|function ob = diag_sp(diag)
%| Construct diagonal "sparse" object, that does D * x = diag(d) * x via d .* x.
%|
%| Copyright 2005-4-6, Jeff Fessler, University of Michigan

if nargin == 1 && streq(diag, 'test'), diag_sp_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if length(varargin) ~= 0
	help(mfilename), error 'bad number of arguments'
end

arg.diag = diag(:);
arg.dim = [1 1] * numel(diag);

%
% build Fatrix object
%
ob = Fatrix(arg.dim, arg, 'forw', @diag_sp_forw, 'back', @diag_sp_back, ...
	'abs', @diag_sp_abs, 'power', @diag_sp_power, ...
	'gram', @diag_sp_gram, 'caller', mfilename);


%
% dot_times_conform()
% (d .* x) operation but reshape d to conform to size of x
%
% in
%	d	[*N 1]
%	x	[(N) (L)]
% out
%	y	[(N) (L)]	repmat(d) .* x
%
function y = dot_times_conform(d, x)
siz = size(x);
if numel(d) == 1 % scalar special case
	y = d * x;
return
end
nd = min(find(numel(d) == cumprod(siz)));
if isempty(nd)
	whos
	fail('dimension mismatch')
end
d = reshape(d, [siz(1:nd) 1]); % [(N)]
d = repmat(d, [ones(1,nd) siz(nd+1:end)]); % [(N) (L)]
y = d .* x;


%
% diag_sp_forw(): y = D * x
%
function y = diag_sp_forw(arg, x)
y = dot_times_conform(arg.diag, x);
%y = repmat(arg.diag, 1, ncol(x)) .* x;


%
% diag_sp_back(): x = D' * y
%
function x = diag_sp_back(arg, y)
x = dot_times_conform(conj(arg.diag), y);
%x = repmat(conj(arg.diag), 1, ncol(y)) .* y;


%
% diag_sp_abs(): |D|
%
function ob = diag_sp_abs(ob, p)
ob.arg.diag = abs(ob.arg.diag);


%
% diag_sp_power(): D .^ p
%
function ob = diag_sp_power(ob, p)
ob.arg.diag = ob.arg.diag .^ p;


%
% diag_sp_gram(): D'*W*D
%
function [ob reuse] = diag_sp_gram(ob, W, reuse, varargin)
if isempty(W)
	ob.arg.diag = abs(ob.arg.diag).^2;
else
	fail('W not implemented yet')
end


%
% diag_sp_test()
%
function diag_sp_test

d = [-2:5]';
d = d + 1i*rand(size(d));
D = diag_sp(d);
Fatrix_test_basic(D, true(size(d)), 'complex', 1)

x = d.^3;
y = D*[x x]; y = y(:,1);
jf_equal(d.*x, y)
x = D'*[y y]; x = x(:,1);
jf_equal(conj(d).*y, x)
y = D.^2*x;
jf_equal(d.^2.*x, y)
y = abs(D)*x;
jf_equal(abs(d).*x, y)

T = build_gram(D, []);
equivs(T * x, D' * (D * x));

d = 7;
D = diag_sp(d);
x = ones(10,1);
jf_equal(d * x, D * x);

