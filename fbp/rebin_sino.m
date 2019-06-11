  function sino2 = rebin_sino(sino1, geom1, geom2, varargin)
%|function sino2 = rebin_sino(sino1, geom1, geom2, varargin)
%|
%| Rebin a sinogram from the geomtry in "geom1" into the geometry in "geom2"
%| both of which were created using sino_geom().
%| The typical use is to convert between fan-beam and parallel-beam.
%|
%| in
%|	sino1	[nb1 na1]	sinogram
%|	geom1	struct		geometry of input sinogram, from sino_geom()
%|	geom2	struct		geometry of ouput sinogram, from sino_geom()
%|
%| option
%|	ob	scalar	1 to create object rather than doing it.  (default: 0)
%|
%| out
%|	sino2	[nb2 na2]	sinogram
%|
%| Copyright 2006-1-18, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sino1, 'test'), rebin_sino_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.ob = 0;
arg = vararg_pair(arg, varargin);

if streq(geom1.type, 'fan') || streq(geom2.type, 'moj')
	sino2 = rebin_fan2par(sino1, geom1, geom2, 'ob', arg.ob);
else
	sino2 = rebin_par2fan(sino1, geom1, geom2, 'ob', arg.ob);
end


% rebin_sino_test()
function rebin_sino_test
rebin_fan2par test
rebin_par2fan test
