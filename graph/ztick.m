 function ztick(arg)
%function ztick(arg)
%	set axis zticks to just end points

if ~nargin
	lim = get(gca, 'zlim');
	if lim(1) == -lim(2)
		lim = [lim(1) 0 lim(2)];
	end
	set(gca, 'ztick', lim)

elseif nargin == 1
	set(gca, 'ztick', arg)
else
	error arg
end
