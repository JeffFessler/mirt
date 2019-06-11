 function [idims_same idims] = fatrix2_same_idims(blocks)
%function [idims_same idims] = fatrix2_same_idims(blocks)
%|
%| check if idims if each block is same;
%| if so, return aggregate idims too as cell

MM = numel(blocks);
b1 = blocks{1};
idims = cell(MM, 1);
idims_same = true;
for mm=1:MM
	bm = blocks{mm};
	if ~isequal(numel(bm.idim), numel(b1.idim))
		idims_same = false;
	end
	idims{mm} = bm.idim;
	idims_same = idims_same && isequal(bm.idim, b1.idim);
end
