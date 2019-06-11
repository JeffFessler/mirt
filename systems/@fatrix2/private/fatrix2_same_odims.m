 function [odims_same odims] = fatrix2_same_odims(blocks)
%function [odims_same odims] = fatrix2_same_odims(blocks)
%|
%| check if odims if each block is same;
%| if so, return aggregate odims too as cell

MM = numel(blocks);
b1 = blocks{1};
odims = cell(MM, 1);
odims_same = true;
for mm=1:MM
	bm = blocks{mm};
	if ~isequal(numel(bm.odim), numel(b1.odim))
		odims_same = false;
	end
	odims{mm} = bm.odim;
	odims_same = odims_same && isequal(bm.odim, b1.odim);
end
