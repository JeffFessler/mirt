 function xtick(arg)
%function xtick(arg)
%|
%| set axis xticks to just end points

if ~nargin
	lim = get(gca, 'xlim');
	if lim(1) == -lim(2)
		lim = [lim(1) 0 lim(2)];
	end
	set(gca, 'xtick', lim)

elseif nargin == 1
	if ischar(arg)
		switch arg
		case 'off'
			set(gca, 'xtick', [])
			set(gca, 'xticklabel', [])
		case 'test'
			xtick_test
		otherwise
			fail 'bug'
		end
	else
		set(gca, 'xtick', sort(arg), 'fontsize', ir_fontsize('tick'))
	end
else
	error arg
end

function xtick_test
im clf, im(eye(7))
%xtick
