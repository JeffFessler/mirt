function varargout = eig(ob, varargin)
% compute eig of object by converting it to full matrix.
% this will be very slow for big objects

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
