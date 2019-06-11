 function ob_out = kronI(M, ob_in)
%function ob_out = kronI(M, ob_in)
% ob_out = kron(eye(M), ob_in)
% useful for making 3D from 2D objects working slice-by-slice

ob_out = fatrix2_kroni(M, ob_in);
