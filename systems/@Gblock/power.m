 function ob = power(ob, sup)
%function ob = power(ob, sup)
%	array "power" method (G.^2) for Gblock class
%	matrix power (mpower) would be for G^2

ob.base = (ob.base) .^ sup;
