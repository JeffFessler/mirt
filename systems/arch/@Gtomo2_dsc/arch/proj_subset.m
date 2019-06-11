 function y = proj_subset(ob, x, ia_start, ia_inc)
%function y = proj_subset(ob, x, ia_start, ia_inc)
%	y = G(:,ia, :,:) * x
%	caution: ia_start is from 1, not 0!

%if ob.is.empty, error empty, end
if ob.apower ~= 1, error 'power not done', end
if ob.is_transpose, error 'transpose not done', end

if ob.is_masked
	x = embed(x, ob.mask);
end

y = wtfmex('dsc,proj', ob.arg', single(x), ob.mask, ...
	int32(ia_start-1), int32(ia_inc), int32(ob.nthread), ob.chat);

y = y(:,ia_start:ia_inc:ob.na);

y = ob.scale * double(y(:));
