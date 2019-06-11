 function out = ir_sparse(in)
%function out = ir_sparse(in)
%|
%| convert matrix to a sparse matrix

if ir_is_octave
	out = sparse(in);
else
	out = sparse(double(in)); % stupid matlab insists on double
end
