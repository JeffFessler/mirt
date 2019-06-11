function varargout = svd(ob, varargin)
% compute svd of object by converting it to full matrix.
% this will be very slow for big objects

ob = full(ob);
[u s v] = svd(ob, varargin{:});

switch nargout
case 1
	varargout = {diag(s)};
case 2
	varargout = {u, s};
case 3
	varargout = {u, s, v};
otherwise
	fail 'eig bug'
end
