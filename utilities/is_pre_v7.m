  function y = is_pre_v7
%|function y = is_pre_v7
%| return 1 if version is before 7, i.e., before release 14.

if ir_is_octave || isfreemat
	y = false;
return
end

[ver ok] = str2num(version('-release'));
if ~ok
	y = false; % 2012a etc
else
	y = ver < 14;
end
