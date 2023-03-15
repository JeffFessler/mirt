 function ob_out = kronI(M, ob_in, varargin)
%|function ob_out = kronI(M, ob_in, varargin)
%|
%| ob_out = kron(eye(M), ob_in)
%|
%| useful for making 3D from 2D objects working slice-by-slice
%|
%| option
%| 'parfor' default 0

arg.parfor = false;
arg = vararg_pair(arg, varargin);

ob_out = fatrix2_kroni(M, ob_in, arg);
