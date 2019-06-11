 function ir_set_gca_fontname(type)
%function ir_set_gca_fontname(type)
%
% ideally the following line would just work:
% set(gca, 'fontname', type)
% but version 7.5 matlab has a bug where
% it sets the font not only of the current axis, but
% also of the legend.  so get all the current font
% names first and then fix them up.

ax = get(gcf, 'children');
for ii=1:length(ax)
	name{ii} = get(ax(ii), 'fontname');
end
set(gca, 'fontname', type)
for ii=1:length(ax)
	if ax(ii) ~= gca
		if ~streq(name{ii}, get(ax(ii), 'fontname'))
			persistent warned
			if isempty(warned)
				warned = true;
				warn 'fixing matlab font bug, hopefully!'
			end
			set(ax(ii), 'fontname', name{ii})
		end
	end
end
