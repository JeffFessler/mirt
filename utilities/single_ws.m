% single_ws
% script to convert all 'double' variables in the workspace to single.
% does not look at elements of structures or in cells

if isvar('ws_ii') || isvar('ws_tmp') || isvar('ws_vname')
	error 'name conflict, cannot do it'
end

ws_tmp = whos;
for ws_ii = 1:numel(ws_tmp)
	if streq(ws_tmp(ws_ii).class, 'double') && ~ws_tmp(ws_ii).sparse
		ws_vname = ws_tmp(ws_ii).name;
		ws_vname = [ws_vname ' = single( ' ws_vname ');'];
		eval(ws_vname)
	end
end
clear ws_ii ws_tmp ws_vname
