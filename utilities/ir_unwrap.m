 function ph = ir_unwrap(ph, varargin)
%function ph = ir_unwrap(ph, varargin)
%|
%| Unwrap phase along all dimensions.
%|
%| 2014-08-19, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(ph, 'test'), ir_unwrap_test, return, end

for id = 1:ndims(ph)
	ph = unwrap(ph, [], id);
end


% ir_unwrap_test
function ir_unwrap_test
p0 = linspace(-pi,2*pi,98)' * linspace(1,2,99);
p0 = p0.'; % octave worked better this way!?
x = exp(1i * p0);
ph = angle(x);
pu = ir_unwrap(ph);
clim = [-1 2]*pi;
im plc 2 2
im(1, p0, 'true', clim), cbar
im(2, ph, 'angle', clim), cbar
im(3, pu, 'unwrapped', clim), cbar
im(4, angle(exp(1i*(pu -p0)))), cbar
