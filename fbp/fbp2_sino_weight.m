  function sino = fbp2_sino_weight(sg, sino, varargin)
%|function sino = fbp2_sino_weight(sg, sino, varargin)
%|
%| Apply sinogram weighting for first step of 2D fan-beam FBP.
%| This matlab version is the backup alternative for users lacking mex routine.
%|
%| in
%|	sg			sino_geom()
%|	ig			image_geom()
%|	sino	[nb,na,(L)]	fan-beam sinogram(s) (line integrals)
%| out
%|	sino	[nb,na,(L)]	weighted sinogram
%|
%| Copyright 2006-4-19 by Jeff Fessler, University of Michigan

if nargin == 1 && streq(sg, 'test'), fbp2_sino_weight_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.chat = 0;
arg = vararg_pair(arg, varargin);

idim = size(sino);
sino = reshape(sino, idim(1), idim(2), []); % [nb,na,*L]
sino = fbp2_sino_weight_do(sino, sg);
sino = reshape(sino, idim); % [nb,na,(L)]


%
% fbp2_sino_weight_do()
%
function sino = fbp2_sino_weight_do(sino, sg)

if isinf(sg.dfs)
	sino = fbp2_sino_weight_flat(sino, ...
		sg.s, sg.dsd, sg.dso, sg.source_offset);
elseif sg.dfs == 0
	sino = fbp2_sino_weight_arc(sino, ...
		sg.s, sg.dsd, sg.dso, sg.source_offset);
else
	error 'only flat and arc done'
end


%
% fbp2_sino_weight_arc()
%
function sino = fbp2_sino_weight_arc(sino, ss, dsd, dso, source_offset);
na = size(sino,2);
nz = size(sino,3);
gam = ss / dsd;
w1 = abs(dso * cos(gam) - source_offset * sin(gam)) / dsd; % 1D weighting
sino = sino .* repmat(w1, [1 na nz]);


%
% fbp2_sino_weight_flat()
%
function sino = fbp2_sino_weight_flat(sino, ss, dsd, dso, source_offset);
na = size(sino,2);
nz = size(sino,3);
gam = atan(ss / dsd);
w1 = abs(dso * cos(gam) - source_offset * sin(gam)) / dsd; % 1D weighting
sino = sino .* repmat(w1, [1 na nz]);


%
% fbp2_sino_weight_test()
%
function fbp2_sino_weight_test
sg = sino_geom('ge1', 'down', 4);
sino = sg.ones;
s1 = fbp2_sino_weight(sg, sino);

im pl 1 2
im(1, s1), cbar
%max_percent_diff(s1,s2)
