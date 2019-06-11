 function ob = ctranspose(ob)
%function ob = ctranspose(ob)
% "ctranspose" method for Fatrix class

ob.is_transpose = ~ob.is_transpose;

% [ob.index1 ob.index2] = deal(ob.index2, ob.index1);
ob.dim = fliplr(ob.dim);
