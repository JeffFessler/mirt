 function ytick(arg)
%function ytick(arg)
%| set axis yticks to just end points, or empty, or specified values

if ~nargin
	lim = get(gca, 'ylim');
	if lim(1) == -lim(2)
		lim = [lim(1) 0 lim(2)];
	end
	set(gca, 'ytick', lim)

elseif nargin == 1
	if ischar(arg) && streq(arg, 'off')
		set(gca, 'ytick', [])
	else   
		set(gca, 'ytick', sort(arg), 'fontsize', ir_fontsize('tick'))
	end

else
	error arg
end
