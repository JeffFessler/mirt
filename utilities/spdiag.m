 function b = spdiag(a, varargin)
%function b = spdiag(a, options)
%| create a sparse matrix with diagonal given by a
%| option:
%|	'nowarn'	do not warn about diag_sp
%| caution: it may be faster to use newer diag_sp() object instead.

if nargin < 1, ir_usage, end

persistent warned
if (~isvar('warned') || isempty(warned)) && nargin < 2
	warning 'spdiag is superceded by diag_sp'
	warned = 1;
end

a = a(:);
a = double(a); % trick: needed because matlab7 doesn't handle single sparse well
n = length(a);
b = sparse(1:n, 1:n, a, n, n, n);
