 function [T, reuse] = build_gram(ob, W, reuse, varargin)
%function [T, reuse] = build_gram(ob, W, reuse, varargin)
% build "gram matrix object" T = G' W G
% in
%	ob		the system matrix G
%	W		typically diag(wi)
%	reuse		stuff returned by previous call that can be reused
%			the next call to save recomputing things.
%			pass empty matrix [] the first time.
%	varargin	options to be passed to object's internal gram maker
% out
%	T		fatrix object for G' W G
%	reuse		output stuff that could be passed back into this
%			routine the *next* time user calls it to save work.

if ~isvar('W'), W = []; end
if ~isvar('reuse'), reuse = []; end

if isempty(ob.handle_gram)
	if isempty(W)
		T = ob' * ob;
	else
		T = ob' * W * ob;
	end
else
	[T, reuse] = feval(ob.handle_gram, ob, W, reuse, varargin{:});
end
