function out = oselect(ob, y)
% recall that imask/omask and idim/odim are swapped on ctranspose
out = fatrix2_select(ob.omask, ob.odim, y);
