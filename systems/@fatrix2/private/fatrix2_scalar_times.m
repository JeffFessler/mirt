 function ob = fatrix2_scalar_times(scalar, ob)
%function ob = fatrix2_scalar_times(scalar, ob)
%|
%| Construct fatrix2 object for: ob = scalar * ob
%|
%| in
%|	scalar	[1]		numeric scalar
%|	ob	fatrix2
%|
%| out
%|	ob	fatrix2		scalar * ob
%|
%| Copyright 2010-12-04, Jeff Fessler, University of Michigan

%if nargin ~= 2, help(mfilename), error(mfilename), end

ob.scale = ob.scale * scalar;
