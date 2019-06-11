 function x = back_subset(ob, y, ia_start, ia_inc)
%function x = back_subset(ob, y, ia_start, ia_inc)
%	x = G(:,ia, :,:)' * y
%	caution: ia_start is from 1, not 0!

if ob.apower ~= 1, error 'power not done', end
if ~ob.is_transpose, error 'transpose needed', end

ia = ia_start:ia_inc:ob.na;
tmp = zeros(ob.nb,ob.na);
tmp(:,ia) = reshape(y, ob.nb, length(ia));
y = tmp;
x = wtfmex('dsc,back', ob.arg', single(y), ob.mask, ...
	int32(ia_start-1), int32(ia_inc), int32(ob.nthread), ob.chat);
if ob.is_masked
	x = x(ob.mask);
end

x = ob.scale * double(x(:));
