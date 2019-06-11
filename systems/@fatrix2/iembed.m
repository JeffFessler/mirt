function out = iembed(ob, x)
% recall that imask/omask and idim/odim are swapped on ctranspose
out = fatrix2_embed(ob.imask, ob.idim, x);
