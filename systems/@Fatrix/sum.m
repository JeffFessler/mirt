 function s = sum(ob)
%function s = sum(ob)
% "sum" method for this class
% sum(M) = 1' * M = (M' * 1)'

s = (ob' * ones(ob.dim(1),1))';
