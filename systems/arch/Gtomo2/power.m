 function ob = power(ob, sup)
%function ob = power(ob, sup)
%	array "power" method (G.^2) for Gtomo class
%	matrix power (mpower) would be for G^2

ob.apower = ob.apower * sup;
ob.scale = ob.scale^sup;
