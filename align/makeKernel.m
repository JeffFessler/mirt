function [kernel kernelg] = makeKernel(deg, m)
%function [kernel kernelg] = makeKernel(deg, m)
%
% Create a filter kernel
%
% in:
%       'deg'           B-spline degree 
%       'm'             magnification factor (even)
%
% out:
%       kernel          FIR interpolation filter coefficients 
%       kernelg         FIR derivative filter coefficients 
%
% Copyright August 2006, Se Young Chun and Jeff Fessler, University of Michigan
switch deg
	case 3
		kernel = [1/6 4/6 1/6];
		kernelg = [0.5 0 -0.5];
		if m > 1
			o = ones([1, m]);
			for i = 1 : deg+1
				kernel = conv(kernel, o);
				kernelg = conv(kernelg, o);
			end
			kernel = kernel / m^deg;
			kernelg = kernelg / m^deg;
		end

	otherwise

end
