 function [C] = CsparsePz(dim)
%function [C] = CsparsePz(dim)
%
% create 3D z-direction penalty for topology-preserving image registration
% assume zero end condition
% in:
%	dim	[n1 n2 n3]	if logical, then just support (aka mask)
% out:
%	C [n1*(n2+1)*n3, n1*n2*n3]  most rows have a 1 and a -1
%
% Copyright 2008-02-07 Se Young Chun, University of Michigan

switch length(dim)
  case 3
	n1 = dim(1);
	n2 = dim(2);
	n3 = dim(3);

	idx = reshape( 1:n1*n2*n3, [n1 n2 n3] ); 
	idx1 = idx(:,:,1:n3-1);
	idx2 = idx(:,:,2:n3);

	C = sparse(1:n1*n2*(n3-1), idx2(:), 1, n1*n2*(n3-1), n1*n2*n3) - ...
		sparse(1:n1*n2*(n3-1), idx1(:), 1, n1*n2*(n3-1), n1*n2*n3);

  otherwise
	printf('error - this dim can not be supported!')

end


