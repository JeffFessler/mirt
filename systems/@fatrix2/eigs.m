function out = eigs(ob, varargin)
% compute largest magnitude eigenvalue of object by calling eigs,
% which in turn probably uses the power method.
% this may be slow for big objects because it is iterative.

warn 'todo: i could not get this to work'
% A = Gdft('mask', true(8,6));
% eigs(A, 1)


if numel(varargin) ~= 1
	fail 'only eigs(A,1) supported'
end
if ~isequal(varargin{1}, 1)
	fail 'only eigs(A,1) supported'
end

fun = @(x) eigs_fatrix2(ob, x);
out = eigs(fun, 1);

%{
ob = full(ob);
[v d] = eig(ob, varargin{:});

switch nargout
case 1
	varargout = {diag(d)};
case 2
	varargout = {v, d};
otherwise
	fail 'eig bug'
end
%}

function out = eigs_fatrix2(ob, x)
pr numel(x)
if numel(x) == 1 % matlab first calls it with scalar "1"
	x = ones(size(ob,2), 1);
end
out = ob * x;
pr size(out)
