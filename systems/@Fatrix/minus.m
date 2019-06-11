 function c = minus(a, b)
%function c = minus(a, b)
% "minus" method for this class

if isa(a, 'Fatrix')
 a = a * eye(a.dim(2));
end
if isa(b, 'Fatrix')
 b = b * eye(b.dim(2));
end
c = a - b;
