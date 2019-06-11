 function s = sum(ob)
%function s = sum(ob)
%	"sum" method for Gtomo2 class
%	sum(G)	= 1' * G	= (G' * 1_nd)'
%	sum(G')	= 1' * G'	= (G * 1_np)'

if ~isempty(ob.index1) | ~isempty(ob.index2)
	error 'sum after subscript reference not done'
end

s = (ob' * ones(ob.dims(1),1))';
