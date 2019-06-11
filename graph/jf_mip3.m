  function out = jf_mip3(vol, varargin)
%|function out = jf_mip3(vol, varargin)
%|
%| create a MIP (maximum intensity projection) mosaic from a 3D volume
%| in	
%|	vol	[]	[nx ny nz] 3d	
%| option
%|	'show'	bool	default: 1 if no output, 0 else
%|	'type'	char	'mip' (default) | 'sum' | 'mid'
%| out
%|	out	[]	[nx+nz ny+nz] 2d	
%|
%| if no output, then display it using im()
%|
%| Copyright 2010-02-22, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(vol, 'test'), jf_mip3_test, return, end

arg.show = ~nargout;
arg.type = 'mip';
arg = vararg_pair(arg, varargin);

if ndims(vol) > 3, fail '3D only', end

switch arg.type
case 'mip'
	fun = @(v,d) jf_mip3_max(v, d);
case 'sum'
	fun = @(v,d) sum(v, d);
case 'mid'
	pn = jf_protected_names;
	fun = @(v,d) pn.mid3(v, d);
otherwise
	fail('unknown type %s', arg.type)
end

xy = fun(vol, 3); % [nx ny]
xz = fun(vol, 2); % [nx nz]
yz = fun(vol, 1); % [ny nz]

xz = permute(xz, [1 3 2]);
zy = permute(yz, [2 3 1])';
nz = size(vol,3);

out = [	xy, xz;
	zy, zeros(nz,nz)];

if arg.show
	im(out)
end

if ~nargout
	clear out
end

function out = jf_mip3_max(in, d)
if isreal(in)
	out = max(in, [], d);
else
	out = max(abs(in), [], d);
end


function jf_mip3_test
ig = image_geom('nx', 64, 'ny', 60, 'nz', 32, 'dx', 1);
vol = ellipsoid_im(ig, [0 0 0 ig.nx/3 ig.ny/4 ig.nz/3 20 0 1], 'oversample', 2);
mip = jf_mip3(vol);
%im(vol), prompt
im(mip), prompt
jf_mip3(vol, 'type', 'sum')
jf_mip3(vol, 'type', 'mid')
