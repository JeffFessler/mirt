  function [T, reuse] = build_gram(ob, W, reuse, varargin)
%|function [T, reuse] = build_gram(ob, W, reuse, varargin)
%| build "gram matrix object" T = A' W A
%| in
%|	ob		the system matrix A
%|	W		typically diag(wi), having size [1 1]*size(ob,1)
%|	reuse		stuff returned by previous call that can be reused
%|			the next call to save recomputing things.
%|			pass empty matrix [] the first time.
%|	varargin	options to be passed to object's internal gram maker
%|
%| out
%|	T		fatrix2 object for A' W A
%|	reuse		output stuff that could be passed back into this
%|			routine the *next* time user calls it to save work.

if ~isvar('W'), W = []; end
if ~isvar('reuse'), reuse = []; end

if ~isempty(ob.handle_gram)
	[T, reuse] = ob.handle_gram(ob, W, reuse, varargin{:});
return
end

if isempty(W)
	T = ob' * ob;
return
end

if ~isa(W, 'fatrix2') % force W to be a fatrix2
	W = Gmatrix(W, 'idim', ob.odim, 'imask', ob.omask, ...
		'odim', ob.odim, 'omask', ob.omask);
end

T = ob' * (W * ob);

%T = fatrix2_build_gram_mixed(ob, W);


%{
% todo: clear all this?


% fatrix2_build_gram_mixed()
% tries to handle array-mode well when ob is fatrix2 and W is something else
function T = fatrix2_build_gram_mixed(ob, W)

arg.ob = ob;
arg.W = W;
arg.ob_omask = ob.omask;

if size(W,2) ~= size(ob,1)
	fail('size(W,2)=%d but size(ob,1)=%d', size(W,2), size(ob,1))
end

T = fatrix2('arg', arg, ...
	'idim', ob.idim, 'imask', ob.imask, ...
	'odim', ob.idim, 'omask', ob.imask, ...
	'forw', @fatrix2_gram_smart_forw);


% fatrix2_gram_smart_forw()
% x [(N)]
% z [(N)]
function z = fatrix2_gram_smart_forw(arg, x)

y = arg.ob.handle_forw(arg.ob.arg, x);
if ~isempty(arg.W)
	if ~isempty(arg.ob_omask)
		y = y(arg.ob_omask);
		y = arg.W * y;
		y = embed(y, arg.ob_omask);
	else
		y = arg.W * y;
	end
end
z = arg.ob.handle_back(arg.ob.arg, y); % [(N)]

%}
