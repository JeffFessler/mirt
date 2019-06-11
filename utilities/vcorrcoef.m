 function c = vcorrcoef(u, v)
%function c = vcorrcoef(u, v)
% correlation coefficient between two vectors

normu = norm(u(:));
normv = norm(v(:));
if normu == 0 && normv == 0
 c = 1;
elseif normu == 0 || normv == 0
 c = nan;
else
 c = u(:)' * v(:) / normu / normv;
end
