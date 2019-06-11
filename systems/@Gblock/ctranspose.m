 function ob = ctranspose(ob)
%function ob = ctranspose(ob)
%	"ctranspose" method for Gtomo2 class

%	transpose base object
base = ob.base;
ob.base = base';
