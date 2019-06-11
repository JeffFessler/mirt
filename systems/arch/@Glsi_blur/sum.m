 function s = sum(ob)
%function s = sum(ob)
%	s = G * 1	or s = G' * 1

if ob.is.empty
	error empty
end

if ob.power ~= 1, error 'sum: power not done', end

if ob.is.subref, error 'sum: subref not done', end

%
%	s = sum(G') = 1' * G' = (G*1)'
%
if ob.is.transpose
	x = double(ob.mask);
	if ob.is.masked
		x = x(ob.mask);
	end
	s = (ob' * x)';

%
%	s = sum(G) = 1' * G = (G' * 1)'
%
else
	y = ones(ob.nb, ob.na);
	s = (ob' * y(:))';
end
