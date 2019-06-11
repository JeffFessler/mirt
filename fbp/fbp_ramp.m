 function [h, nn] = fbp_ramp(type, n, ds, dsd)
%function [h, nn] = fbp_ramp(type, n, ds, dsd)
%|
%| 'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
%| This sampled band-limited approach avoids the aliasing that would be
%| caused by sampling the ramp directly in the frequency domain.
%|
%| in
%|	type		'arc' (3rd generation CT) or 'flat' (for parallel too)
%|	n	(int)	# of samples (must be even)
%|	ds	(real)	sample spacing (in distance units, e.g., cm)
%|	dsd	(real)	source-to-detector distance, for 'arc' case
%| out
%|	h	[n]	samples of band-limited ramp filter
%|
%| Copyright 2005-6-8, Jeff Fessler, University of Michigan

if nargin == 1 && streq(type, 'test'), fbp_ramp_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

if 2 * floor(n/2) ~= n, error 'n must be even', end

switch type
case 'arc'
	if nargin ~= 4, help(mfilename), error 'need 4 args', end
	[h nn] = ramp_arc(n, ds, dsd);
case 'flat'
	if nargin ~= 3 && ~isempty(dsd), warn('only need 3 args for flat'), end
	[h nn] = ramp_flat(n, ds);
otherwise
	error 'bad fan type'
end


% ramp_arc()
% 'ramp-like' filter for arc fan beam, p.83 of Kak
%
function [h, nn] = ramp_arc(n, ds, dsd)
if n/2 * ds / dsd > 0.9 * pi/2
	printm('angle is %g degrees: too large', rad2deg(n/2 * ds / dsd))
	error 'physically impossible arc geometry'
end
nn = [-(n/2):(n/2-1)]';
h = zeros(size(nn));
h(nn==0) = 1 / (4 * ds^2);
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * dsd * sin(nn(odd) * ds / dsd)).^2;


% ramp_flat()
% 'ramp-like' filter for flat fan beam
% (usual ramp filter for parallel-beam case)
%
function [h, nn] = ramp_flat(n, ds)
nn = [-(n/2):(n/2-1)]';
h = zeros(size(nn));
h(n/2+1) = 1 / 4; 
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * nn(odd)).^2;
h = h / ds^2;


% fbp_ramp_test
function fbp_ramp_test
nb = 256;
ds = 1;
dsd = 100;
[h1 nn] = fbp_ramp('arc', nb, ds, dsd);
[h2 nn] = fbp_ramp('flat', nb, ds);
max_percent_diff(h1, h2)
if im
	plot(nn, h1, 'o', nn, h2, '+')
	xlabel 'n', ylabel 'h[n]'
	axisx(-10, 10)
	legend('arc', 'flat')
end
