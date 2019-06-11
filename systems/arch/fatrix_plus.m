 function ob = fatrix_plus(varargin)
%function ob = fatrix_plus(fatrix_objects)
%
% Construct Fatrix object consisting of sum of Fatrix objects, i.e.,
% T1 + T2 + ... + Tn
%
% in
%	fatrix_objects	one or more fatrix objects
%
% out
%	ob		T1 + ... + Tn
%
% Copyright 2007-10-16, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), fatrix_plus_test, return, end

arg.blocks = varargin;

% options
%arg = vararg_pair(arg, varargin);

arg.dim = size(arg.blocks{1});
for ii=2:length(arg.blocks)
	if ~isequal(size(arg.blocks{2}), arg.dim)
		error 'size mismatch'
	end
end

% build Fatrix object
ob = Fatrix(arg.dim, arg, 'caller', 'fatrix_plus', ...
	'forw', @fatrix_plus_forw, 'back', @fatrix_plus_back);


%
% fatrix_plus_forw(): y = G * x
%
function y = fatrix_plus_forw(arg, x)

y = arg.blocks{1} * x;
for ii=2:length(arg.blocks)
	tmp = arg.blocks{ii} * x;
	y = y + tmp;
end


%
% fatrix_plus_back(): x = G' * y
%
function x = fatrix_plus_back(arg, y)

x = arg.blocks{1}' * y;
for ii=2:length(arg.blocks)
	tmp = arg.blocks{ii}' * y;
	x = x + tmp;
end


% fatrix_plus_test()
function fatrix_plus_test

A1 = diag_sp(1:3);
A2 = diag_sp(4:6);
x = [7:9]';
A = fatrix_plus(A1, A2);
y0 = A * x;
ys = A1 * x + A2 * x;
jf_equal(y0, ys)

A = block_fatrix({A1, A2}, 'type', 'col');
T = build_gram(A, []);
y0 = A * x;
ys = [A1 * x; A2 * x];
jf_equal(y0, ys)
