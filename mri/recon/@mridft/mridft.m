function ob = mridft(kx, ky, t, we, xval, yval, fov, N)
%function ob = mridft(kx, ky, t, we, xval, yval, fov, N)
%	Construct MRI object, which can do Ax and A'y operations
%       Brad Sutton    University of Michigan


%	default object
ob.kx = 0;	% should these be []'s
ob.ky = 0;
ob.t = 0;
ob.we = 0;
ob.xval = 0;
ob.yval = 0;
ob.is.empty	= true;
ob.is.transpose = false;
ob.fov = 0;  %FOV in cm's
ob.n = 0;
%ob.version = 1.0;

if nargin == 0
	ob = class(ob, 'mridft');
	return
end

if isa(kx, 'mridft')
	ob = kx;
	return
end

if nargin ~= 8
	help mridft
	error nargin
end

	%	fill object
	ob.kx = kx;
	ob.ky = ky;
	ob.t = t;
	ob.we = we;
	ob.xval = xval;
	ob.yval = yval;
        ob.fov = fov;
        ob.n = N;        

	ob.is.empty	= false;

%	ob.m = size(we);	% image size
%	ob.n = size(we);

	ob = class(ob, 'mridft');






