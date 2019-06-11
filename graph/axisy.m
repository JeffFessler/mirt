 function axisy(arg1, arg2)
%function axisy(vals)
%function axisy(val1, val2)
% set y range of plot axis

if nargin < 1, help axisy, error args, end
if nargin == 1, vals = arg1; end
if nargin == 2, vals = [arg1 arg2]; end

if ischar(arg1) && streq(arg1, 'tight')
	xlims = get(gca, 'xlim');
	axis tight
	xlim(xlims)
return
end

ylim(vals)
