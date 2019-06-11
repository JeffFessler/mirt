 function ob = embed_fatrix(mask)
%function ob = embed_fatrix(mask)
% Under Construction: perhaps not needed?
% in
%	mask	[dim]	logical support mask
% out
%	ob	[*dim,np] np = sum(mask(:))
arg.nn = numel(mask);
arg.np = sum(mask(:));
dim = [arg.nn arg.np];
ob = fatrix(dim, 'forw', @embed_fatrix_forw, 'back', @embed_fatrix_back);

function y = embed_fatrix_forw(arg, x)
todo

function x = embed_fatrix_back(arg, y)
todo

y = embed_mult(fun, arg, is_transpose, x, 1, 1, ...
 mask, np, odim, is_array, varargin)

