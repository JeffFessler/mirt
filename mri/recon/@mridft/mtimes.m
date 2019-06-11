 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%       Brad Sutton, University of Michigan

if a.is.empty
	error empty
end



if ~a.is.transpose
	vo = mvmult(a.kx, a.ky, a.t, a.we, vi(:), a.xval, a.yval, 0, a.fov, a.n);
else
	vo = mvmult(a.kx, a.ky, a.t, a.we, vi(:), a.xval, a.yval, 1, a.fov, a.n);
end
