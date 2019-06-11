 function axisx(arg1, arg2)
%function axisx(vals) - redundant with xlim()
%function axisx(val1, val2)
%|
%| set x range of plot axis

if nargin < 1, help axisx, error args, end
if nargin == 1, vals = arg1; end
if nargin == 2, vals = [arg1 arg2]; end

if ischar(arg1) && streq(arg1, 'tight')
	ylims = get(gca, 'ylim');
	axis tight
	ylim(ylims)
return
end

xlim(vals)
