  function [sino z_orbit_start] = helix_rebin(cg, ig, proj)
%|function [sino z_orbit_start] = helix_rebin(cg, ig, proj)
%|
%| interface to helix_rebin_mex mex file for rebinning helical CT data
%| into one fan-beam sinogram per image slice
%|
%| in
%|	cg
%|	ig
%|	proj	[ns nt na]
%|
%| out
%|	sino	[ns na_rebin nz]
%|	z_orbit_start [nz]		orbit_start for each sinogram
%|
%| 2010-07-28, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

nthread = 1; % for now only 1 thread supported

i = @(x) int32(x);
s = @(x) single(x);

[sino z_orbit_start] = helix_rebin_mex('helix,rebin,ssrb', ...
	i(ig.nz), s(ig.dz), s(ig.offset_z), ...
	s(cg.dso), s(cg.dsd), s(cg.dfs), ...
	s([cg.ds cg.dt]), s([cg.offset_s cg.offset_t]), ...
	s(cg.pitch), s(cg.source_z0), ...
	s(cg.orbit), s(cg.orbit_start), ...
	s(proj), i(nthread));
