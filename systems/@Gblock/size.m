 function dim = size(ob, varargin)
%function dim = size(ob, varargin)
% "size" method for Gblock class

base = ob.base;

dim = size(base, varargin{:});

if ob.i_block
	if length(varargin), error 'not done', end

	if rem(dim(1), ob.nblock)
		printf('full dim(1)=%d nblock=%d', dim(1), ob.nblock)
		error('fix: cannot report size in non-multiple case!')
	end
	dim(1) = dim(1) / ob.nblock;
end

% old 2d way
%	if isfield(struct(base), 'nb')
%		nb = base.nb;
%		na = base.na;
%	else
%		return
%	end
%	ia = ob.i_block:ob.nblock:na;
%	dim(1) = nb * length(ia);
