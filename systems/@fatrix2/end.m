 function out = end(ob, k, n)
%function out = end(ob, k, n)
%
% handle calls of the following forms, for example:
% ob(:,end)
% ob(end,:)
% ob(end)
%
% see 'doc end' for more explanation

switch n
case 1
	if k ~= 1
		fail 'bug'
	end

	out = prod(size(ob));

case 2
	out = size(ob, k);

otherwise
	fail 'unknown "n"'
end
