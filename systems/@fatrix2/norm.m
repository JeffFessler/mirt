function out = norm(ob, varargin)

if prod(size(ob)) > (2^12)^2 % 64^2 by 64^2
	warn 'norm of large object can be very slow'
end
out = full(ob);
out = norm(out, varargin{:});
